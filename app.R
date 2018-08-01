
library(shiny)
library(ChemmineR)
library(webchem)
library(ChemmineOB)
library(rcdk)
library(caret)
library(gplots)
library(ggplot2)
library(lattice)
library(klaR)
library(latticeExtra)
library(e1071)
library(randomForest)
library(ROCR)
library(party)
library(kernlab)
library(MASS)
library(dplyr)
library(reshape2)
library(reshape)
library(rpart)				    
library(rattle)					
library(rpart.plot)				
library(RColorBrewer)			
library(PerformanceAnalytics)
library(party)		
library(partykit)				
library(tree) 
library(e1071) 
library(xtable)
library(pROC)
library(AppliedPredictiveModeling)
library(randomForest)
if (interactive()) {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        h4(strong('Machine learning models available')),
        # the actioButton called rpart which is the name of the variable you need to use in the server component
        actionButton('rf', label = 'Random Forest', icon("tree-conifer", lib="glyphicon"),
                     style="color: #fff; background-color: #33cc33; border-color: #2e6da4"),
        actionButton('nb', label = 'k-Nearest Neighbors', icon("random", lib="glyphicon"),
                     style="color: #fff; background-color: #0066ff; border-color: #2e6da4"),
        actionButton('svm', label = 'Support Vector Machine ', icon("resize-full", lib="glyphicon"),
                     style="color: #fff; background-color: #ffa500; border-color: #2e6da4"),
        tags$hr(),
        h4(numericInput("ratio", "Specify the % for training sample", value=50/100, min = 50/100, max = 90/100, step=0.1)),
        h4(strong('Select the input format')),
        radioButtons("colour",label=NA,choices=c("InChI","InChI KEY","CAS ID","PUBCHEM SID/CID"),selected="InChI"),
        h4(fileInput("filed", strong("Upload CSV File")))
      ),
      mainPanel(
        tableOutput("contents"),
        tabsetPanel(
          tabPanel("Predictive Modeling", verbatimTextOutput("headi")), 
          #tabPanel("Info on the machine learning models", tableOutput("results")), 
          tabPanel("ROC Curve", plotOutput('ploti')),
          tabPanel("Ensemble Modeling", verbatimTextOutput("ensem"))
          #tabPanel("Info on the machine learning models", verbatimTextOutput("confusion"))
        ),
        strong(helpText("Prediction Results Using Random Forest")),
        verbatimTextOutput("outputs"),
        strong(helpText("Prediction Results Using Support Vector Machine")),
        verbatimTextOutput("outputs2"),
        strong(helpText("Prediction Results Using kNN")),
        verbatimTextOutput("outputs3")
      )
    )
  )
  server <- function(input, output,session) {
    atad <- reactive({
      ControlParamteres <- trainControl(method = "cv",number =5,savePredictions = TRUE,classProbs = TRUE)
      parameterGrid <- expand.grid(mtry=c(2,3,4))
      ctrl <- trainControl(method="repeatedcv",repeats = 3,classProbs=TRUE,summaryFunction = twoClassSummary)
      trControl <- trainControl(method  = "cv",number= 10,classProbs = TRUE,summaryFunction = twoClassSummary)
      if (input$colour == "PUBCHEM SID/CID"){
        filej <- input$filed
        if(is.null(filej)){return()}
        res <-read.table(file=filej$datapath,sep="\t",header = TRUE)
        cas <- res[,2]
        smiles <- NULL
        for (x in 1:length(cas)){
          smiles <- c(smiles, pc_prop(as.character(cas[x]),properties='CanonicalSMILES')[[2]])
        }
        res$SMILES <- smiles
        write.table(res[!is.na(res$SMILES),], "smilesres.csv", sep="\t", row.names = F)
        write.table(res[is.na(res$SMILES),], "random_NA1.csv", sep="\t", row.names=F)
      }
      else  if(input$colour == "InChI"){
        filej <- input$filed
        if(is.null(filej)){return()}
        res <-read.csv(file=filej$datapath,sep="\t",header = TRUE)
        inchin <-res[,4]
        inchin
        extractinchi <-write.csv(inchin, "extractedinchi.csv",row.names=T)
        fromfile <-read.csv("extractedinchi.csv", sep="\t", header=T)
        fromfile
        y <-convertFormatFile("inchi","smiles","extractedinchi.csv","smiles.csv")
        tofile<-read.table("smiles.csv",header = FALSE)
        smiles <- NULL
        smiles <-tofile[,1]
        smiles
        res$SMILES <- NA
        res$SMILES<-smiles
        write.table(res[(res$SMILES),], "smilesres.csv", sep="\t", row.names=F)
      }
      else if (input$colour == "CAS ID" |  input$colour =="InChI KEY"){
        filej <- input$filed
        if(is.null(filej)){return()}
        res <-read.table(file=filej$datapath,sep="\t",header = TRUE)
        cas <- res[,3]
        smiles <- NULL
        for (x in 1:length(cas)){
          smiles <- c(smiles, cir_query(as.character(cas[x]), 'smiles')[[1]])
        }
        res$SMILES <- smiles
        write.table(res[!is.na(res$SMILES),], "smilesres.csv", sep="\t", row.names = F)
        write.table(res[is.na(res$SMILES),], "random_NA1.csv", sep="\t", row.names=F)
      }
      data<-read.csv("smilesres.csv", sep="\t", header=T)
      output$headi <- renderPrint({
        pls <- data[1:5,c(2,5)]
        head(pls)
      })
      #Parse SMILES
      data_filt <- data[!is.na(data$SMILES),]
      mols <- parse.smiles(as.character(data_filt$SMILES))
      names(mols) <-  data_filt$ID
      #Add CDK descriptors
      descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
      CDKdescs <- eval.desc(mols, descNames, verbose=T)
      CDKdescs$ID <- rownames(CDKdescs)
      data_descs <- merge(data_filt, CDKdescs, by="ID")
      write.table(data_descs, "descriptors.csv", sep="\t", row.names=F, quote=F)
      #Data Preprocessing using caret package (This is optional - might not be best to do this if small number of samples per class and sparse PLS methods can handle correlated predictors)
      data.mat <- data_descs
      rownames(data.mat) <- data.mat$ID
      data.mat <- data.mat[,-c(1:5)] #N.B. Change based on col positions
      #Remove zero and near-zero variance predictors
      print(paste("Starting with", ncol(data.mat), "descriptors. Beginning filtering process...", sep=" "))
      nzv <- nearZeroVar(data.mat, saveMetrics= TRUE)
      write.table(nzv, "descriptor_variance_summary.csv", sep="\t", quote=F, col.names=NA)
      nzv <- nearZeroVar(data.mat)
      filt.data.mat <- data.mat[,-nzv]
      print(paste("Removed", length(nzv), "near-zero variance predictors...", ncol(filt.data.mat), "descriptors remaining.", sep=" "))
      #removal of highly correlated predictors (likely repeated descriptors from multiple algorithms)
      descrCor <-  cor(na.omit(filt.data.mat))
      highCor <- findCorrelation(descrCor, cutoff = .999, verbose=T)
      filt.data.mat <- filt.data.mat[,-highCor]
      print(paste("Removed", length(highCor), "highly correlated predictors...", ncol(filt.data.mat), "descriptors remaining.", sep=" "))
      write.table(filt.data.mat, "descriptors_filtered.csv", sep="\t", quote=F, col.names=NA)
      descri<- read.csv("descriptors.csv", sep="\t", header=T)
      descri <- descri[,colSums(is.na(descri))<nrow(descri)]
      r<-as.numeric(input$ratio)
      ind <- sample(2, nrow(descri), replace = TRUE, prob=c(r,1-r))
      trainset = descri[ind==1,]
      testset = descri[ind==2,]
      modelRandom <- train(PUBCHEM_ACTIVITY_OUTCOME~., data = trainset,method ="rf",trControl = ControlParamteres,tuneGrid=parameterGrid)
      predictions<-predict(modelRandom,testset)
      actual <-testset$PUBCHEM_ACTIVITY_OUTCOME
      testing1<-table(predictions=predictions,actual=testset$PUBCHEM_ACTIVITY_OUTCOME)
      output$outputs <- renderPrint({
        testing1
      })
      svm_model <- svm(PUBCHEM_ACTIVITY_OUTCOME~., data=trainset,type="C-classification",kernel="linear", cost=1, cross=10)
      pred <- predict(svm_model,testset)
      # Model accuracy
      testing2<-table(predictions=pred,actual=testset$PUBCHEM_ACTIVITY_OUTCOME)
      output$outputs2 <- renderPrint({
        testing2
      })
      knnFit <- train(PUBCHEM_ACTIVITY_OUTCOME~ ., data = trainset, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 10)
      knnPredict <- predict(knnFit,testset)
      tak<-table(knnPredict,testset$PUBCHEM_ACTIVITY_OUTCOME)
      output$outputs3 <- renderPrint({
        tak
      })
      output$ploti <- renderPlot({
        predictionsrf<-predict(modelRandom,testset,type='prob')
        pred.rf<-prediction(predictionsrf[,2],testset$PUBCHEM_ACTIVITY_OUTCOME)
        rf.auc<-performance(pred.rf,'tpr','fpr')
        predictknn<-predict(knnFit,testset,type='prob')
        pred.knn<-prediction(predictknn[,2],testset$PUBCHEM_ACTIVITY_OUTCOME)
        knn.auc<-performance(pred.knn,'tpr','fpr')
        svmmodel<-svm(PUBCHEM_ACTIVITY_OUTCOME~., data=trainset, method="C-classification",
                      kernel="radial", gamma = 0.01, cost = 100,cross=5, probability=TRUE)
        svmmodel.predict<-predict(svmmodel,subset(testset,select=-PUBCHEM_ACTIVITY_OUTCOME),decision.values=TRUE)
        svmmodel.probs<-attr(svmmodel.predict,"decision.values")
        svmmodel.class<-predict(svmmodel,testset,type="class")
        svmmodel.labels<-testset$PUBCHEM_ACTIVITY_OUTCOME
        #roc analysis for test data
        svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
        svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
        svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
        rfauc<-performance(pred.rf,"auc")@y.values[[1]]
        nbauc<-performance(pred.knn,"auc")@y.values[[1]]
        svmauc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
        Methods<- c('Random Forest','kNN','SVM')
        AUCScore<-c(rfauc,nbauc,svmauc)
        data.frame(Methods,AUCScore)
        #PLot and adding legend 
        plot(rf.auc,col='red',lty=1,main='ROC Curve Comparison of Random Forest V/s kNN V/s SVM')
        plot(knn.auc,col='green',add=TRUE,lty=1)
        plot(svmmodel.performance,col='black',add=TRUE,lty=1)
        L<-list(bquote("Random Forest"== .(rfauc)), bquote("kNN"== .(nbauc)),bquote("SVM"== .(svmauc)))
        legend("bottomright",legend=sapply(L, as.expression),col=c('red','green','black'),lwd=2,bg="gray",pch=14,cex=0.6)
      })
    })
    output$contents <- renderTable({
      if(is.null(atad())){return ()}
      atad()
    })
    
  }
  shinyApp(ui, server)
}
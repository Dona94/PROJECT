library(shiny)
library(ChemmineR)
library(webchem)
library(fscaret)
library(ChemmineOB)
library(rcdk)
library(ensembleR)
library(caretEnsemble)
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
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }"
        ),
        h4(strong('Machine learning models available')),
        # the actioButton called rpart which is the name of the variable you need to use in the server component
        actionButton('svm', label = 'Treebag', icon("leaf", lib="glyphicon"),
                     style="color: #fff; background-color: #ffa500; border-color: #2e6da4"),
        actionButton('rf', label = 'Random Forest', icon("tree-conifer", lib="glyphicon"),
                     style="color: #fff; background-color: #33cc33; border-color: #2e6da4"),
        actionButton('nb', label = 'k-Nearest Neighbors', icon("random", lib="glyphicon"),
                     style="color: #fff; background-color: #0066ff; border-color: #2e6da4"),
        
       
        
        h4(numericInput("ratio", "Specify the % for training set", value=50/100, min = 50/100, max = 90/100, step=0.1)),
        h4(strong('Select the input format')),
        radioButtons("colour",label=NA,choices=c("InChI","InChI KEY","CAS ID","PUBCHEM SID/CID"),selected="InChI"),
        h4(fileInput("filed", strong("Upload CSV File")))
      ),
      mainPanel(
        tableOutput("contents"),
        tabsetPanel(
          #tabPanel("Info on the machine learning models", tableOutput("results")), 
          tabPanel("Predictive Modeling-RF,kNN,treebag", verbatimTextOutput("outputs"),verbatimTextOutput("outputs2"),
                   verbatimTextOutput("outputs3")), 
          tabPanel("Confusion Matrix-RF,kNN,treebag", verbatimTextOutput("crf"),verbatimTextOutput("csvm"),
                   verbatimTextOutput("cknn")), 
          tabPanel("ROC Curve", plotOutput('ploti')),
          tabPanel("Descriptors", plotOutput('des'))
          
          #tabPanel("Info on the machine learning models", verbatimTextOutput("confusion"))
        )
      )
    )
  )
  server <- function(input, output,session) {
    atad <- reactive({
      ControlParamteres <- trainControl(method = "cv",number =5,savePredictions = TRUE,classProbs = TRUE)
      parameterGrid <- expand.grid(mtry=c(2,3,4))
      fitControl <- trainControl(method = "cv",number = 5,savePredictions = 'final',
                                 classProbs = T)
      ctrl <- trainControl(method="repeatedcv",repeats = 3,classProbs=TRUE,summaryFunction = twoClassSummary)
      trControl <- trainControl(method  = "cv",number= 10,classProbs = TRUE,summaryFunction = twoClassSummary)
      outcomeName <- 'PUBCHEM_ACTIVITY_OUTCOME'
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
      descri<- read.csv("descriptors.csv", sep="\t", header=T)
      descri <- subset(descri, select = -c(ID, PUBCHEM_SID,InChI,SMILES))
      descri <- descri[,colSums(is.na(descri))<nrow(descri)]
      predictors <- setdiff(names(descri),outcomeName)
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
      model_lr<-train(trainset[,predictors],trainset[,outcomeName],method='treebag',trControl=fitControl,tuneLength=3)
      #Predicting using knn model
      testset$pred_lr<-predict(object = model_lr,testset[,predictors])
      predtbag<-predict(model_lr,testset)
      testing2 <- table(predtbag,testset$PUBCHEM_ACTIVITY_OUTCOME)
      output$outputs2 <- renderPrint({
        testing2
      })
      knnFit <- train(PUBCHEM_ACTIVITY_OUTCOME~ ., data = trainset, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 10)
      knnPredict <- predict(knnFit,testset)
      tak<-table(knnPredict,testset$PUBCHEM_ACTIVITY_OUTCOME)
      output$outputs3 <- renderPrint({
        tak
      })
      rf <-confusionMatrix(testset$PUBCHEM_ACTIVITY_OUTCOME,predictions)
      output$crf <- renderPrint({
        rf
      })
     tbg <- confusionMatrix(testset$PUBCHEM_ACTIVITY_OUTCOME,testset$pred_lr)
      output$csvm <- renderPrint({
      tbg
      })
      knn <-confusionMatrix(testset$PUBCHEM_ACTIVITY_OUTCOME,knnPredict)
      output$cknn <- renderPrint({
        knn
      })
      output$ploti <- renderPlot({
        predictionsrf<-predict(modelRandom,testset,type='prob')
        predictionsrf
        pred.rf<-prediction(predictionsrf[,2],testset$PUBCHEM_ACTIVITY_OUTCOME)
        pred.rf
        rf.auc<-performance(pred.rf,'tpr','fpr')
        rf.auc
        predictknn<-predict(knnFit,testset,type='prob')
        pred.knn<-prediction(predictknn[,2],testset$PUBCHEM_ACTIVITY_OUTCOME)
        knn.auc<-performance(pred.knn,'tpr','fpr')
        predtbag<-predict(model_lr,testset,type='prob')
        predtbag<-prediction(predtbag[,2],testset$PUBCHEM_ACTIVITY_OUTCOME)
        tbag.auc<-performance(predtbag,'tpr','fpr')
        rfauc<-performance(pred.rf,"auc")@y.values[[1]]
        nbauc<-performance(pred.knn,"auc")@y.values[[1]]
        tbgauc<-performance(predtbag,"auc")@y.values[[1]]
        Methods<- c('Random Forest','kNN','TreeBag')
        AUCScore<-c(rfauc,nbauc,tbgauc)
        data.frame(Methods,AUCScore)
        #PLot and adding legend 
        plot(rf.auc,col='red',lty=1,main='ROC Curve Comparison of RF, SVM, treebag')
        plot(knn.auc,col='green',add=TRUE,lty=1)
        plot(tbag.auc,col='black',add=TRUE,lty=1)
        L<-list(bquote("Random Forest"== .(rfauc)), bquote("Naive Bayes"== .(nbauc)),bquote("Tree Bag"== .(tbgauc)))
        legend("bottomright",legend=sapply(L, as.expression),col=c('red','green','black'),lwd=2,bg="gray",pch=20,cex=1)
      })
      desc<- read.csv("descriptors.csv", sep="\t", header=T)
      desc <- subset(desc, select = -c(ID, PUBCHEM_SID,InChI,SMILES))
      desc <- desc[,colSums(is.na(desc))<nrow(desc)]
      #desc$PUBCHEM_ACTIVITY_OUTCOME <- factor(ifelse(desc$PUBCHEM_ACTIVITY_OUTCOME=='Active', 1,0))
      fs<-desc %>% select(-PUBCHEM_ACTIVITY_OUTCOME,everything())
      tDummy <- dummyVars("~.",data=fs, fullRank=F)
      tDF <- as.data.frame(predict(tDummy,fs))
      ind<-sample(2,nrow(fs),replace=TRUE,prob=c(0.8,0.2))
      trainset<-fs[ind==1,]
      testset<-fs[ind==2,]
      fsModels <- c("rf","knn","treebag") 
      myFS<-fscaret(trainset, testset, preprocessData=TRUE,
                    Used.funcRegPred = fsModels, with.labels=TRUE,
                    supress.output=FALSE,no.cores=1)
      results <- myFS$VarImp$matrixVarImp.MSE
      results$Input_no <- as.numeric(results$Input_no)
      results <- results[c("SUM","SUM%","ImpGrad","Input_no")]
      myFS$PPlabels$Input_no <- as.numeric(rownames(myFS$PPlabels))
      results <- merge(x=results, y=myFS$PPlabels, by="Input_no", all.x=T)
      results <- results[c('Labels', 'SUM')]
      results <- subset(results,results$SUM !=0)
      results <- results[order(-results$SUM),]
      print(results)
      results1 = head(results, 10)
      output$des <- renderPlot({
        p <- ggplot(results1, aes(x=Labels, y=SUM)) +
          geom_bar(stat="identity", fill="#53cfff") +
          coord_flip() +
          theme_light(base_size=16) +
          xlab("Predictors") +
          ylab("") +
          ggtitle("Top ten descriptors predictive of activity") +
          theme(plot.title=element_text(size=12))
        p
      })
    })
    output$contents <- renderTable({
      if(is.null(atad())){return ()}
      atad()
    })
    
  }
  shinyApp(ui, server)
}





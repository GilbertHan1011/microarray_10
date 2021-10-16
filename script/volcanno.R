library(ggrepel)
plotData <-  annoTab%>%
  mutate(logP=-log10(adj.P.Val))%>%
  mutate(type = ifelse((logP>1.3&logFC > 1), "UP",
                       ifelse((logP>1.3&logFC < -1),"DOWN","NONE")))%>%filter(logP>0)%>%
  select(logP, logFC, type, SYMBOL,PROBEID)
prePlot <- data%>%filter(logP>0)
g <- ggplot(plotData, aes(x = logFC, y = logP, color = type,label=SYMBOL))+ 
  geom_point(size = 1.3)+
  scale_color_manual(values=c('#0000FF', 'grey','red'))+ 
  theme_classic()+scale_x_continuous(name="log2FoldChange")+
  scale_y_continuous(name="-log10Pvalue") + 
  labs(title = "treatment_vs_control")+
  theme(plot.title = element_text(hjust=0.5,size=16, face="bold"))+
  geom_vline(xintercept = 1, linetype="dashed")+ geom_vline(xintercept = -1, linetype="dashed")+
  geom_hline(yintercept = 1.3, linetype="dashed")+
  theme(plot.title = element_text(hjust=0.5,size=16, face="bold"))+
  theme(text = element_text(size= 20))+
  theme(axis.line = element_line(size=1, colour = "black"))+
  theme(axis.text.x = element_text(size = 15, family = "myFont",face = "bold"), axis.text.y = element_text(size = 15, family = "myFont",face = "bold"))+
  theme(legend.position="none")
g
labelGene <- c("Mt1")
g+geom_text_repel(data = subset(plotData, SYMBOL%in%labelGene & type!="NONE"),box.padding=2,size=6,color="black")

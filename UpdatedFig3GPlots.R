ggplot(mPerctMeltedAnnot,aes(x = as.numeric(FA_Score),y=value))+
  geom_point(aes(col=Recov_LC))+
  theme_classic()+geom_smooth(method = "lm")+xlab("FAS-Score")+ylab("% of cells in CD14 Monocytes cluster per donor")+ #
  facet_wrap(~Var2,scales = "free")+stat_cor(method = "spearman")+scale_color_manual(values = c("acute" = "pink","mild/moderate" = "firebrick2","Recovered"="cyan","severe"="darkred"))+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",
        strip.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(colour = "black",fill = "transparent",linetype = "solid", size = 0.5),
        panel.border = element_rect(colour = "black",fill = "transparent",linetype = "solid", size = 0.5))

ggplot(mPerctMeltedAnnot,aes(x = as.numeric(pO2),y=value))+
  geom_point(aes(col=Recov_LC))+
  theme_classic()+geom_smooth(method = "lm")+xlab("pO2")+ylab("% of cells in CD14 Monocytes cluster per donor")+
  facet_wrap(~Var2,scales = "free")+stat_cor(method = "spearman")+scale_color_manual(values = c("acute" = "pink","mild/moderate" = "firebrick2","Recovered"="cyan","severe"="darkred"))+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",
        strip.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(colour = "black",fill = "transparent",linetype = "solid", size = 0.5),
        panel.border = element_rect(colour = "black",fill = "transparent",linetype = "solid", size = 0.5))

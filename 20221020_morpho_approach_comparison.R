d_old <- as.data.frame(morphology_old_noNA)
ggplot(d_old, aes(x=body_length, y=pVII_length, label=rownames(d_old)))+
  geom_point() +
  geom_text(size=2)+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("old data")

d <- as.data.frame(morphology_noNA)
ggplot(d, aes(x=body_length, y=pVII_length, label=rownames(d)))+
  geom_point() +
  geom_text(size=2)+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("imputed")

d_top3 <- as.data.frame(morpho_AT3)
ggplot(d_top3, aes(x=body_length, y=pVII_length, label=rownames(d_top3)))+
  geom_point() +
  geom_text(size=2)+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("top3only")
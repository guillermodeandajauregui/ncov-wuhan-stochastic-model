plot.model(data = mis_modelos$`1`, log = "", title = "prueba1")
plot.model(data = mis_modelos[1], log = "", title = "prueba2")
plot.model(data = mis_modelos[[5]], log = "", title = "prueba3")

lapply(X = 1:5, FUN = function(i){
  plot.model(data = mis_modelos[[i]], log = "", title = escenarios$Description[[i]])
})

plot.ensemble(mis_modelos)

get.range(mis_modelos)

mis_modelos$`1`

mis_modelos[[1]] %>% bind_rows() %>% select("C") %>% max()

mis_modelos %>% 
  lapply(FUN = function(i){
    mini <- i %>% bind_rows() %>% select("C") %>% min()
    maxi <- i %>% bind_rows() %>% select("C") %>% max()
    return(data.frame(minimum = mini,
                      maximum = maxi))
  }) %>% bind_rows(.id = "escenario")

mis_out <- lapply(mis_modelos, get.range)
names(mis_out) <- paste0("escenario_", 1:5)

plot.ensemble(x = mis_out, plausible = 1)

mis_out %>% bind_rows() %>% cbind(data.frame(parameter = c("min", "max"))) %>% pivot_longer(cols = -parameter) %>% 
  rename(escenario = name, casos_esperados = value) %>% 
  mutate(medidas = ifelse(escenario%in%c("escenario_1", "escenario_3"), "no", "si")) %>% 
  ggplot(aes(x = escenario, y = casos_esperados, colour = medidas)) + 
  geom_line() + 
  theme_minimal() + 
  scale_colour_manual(values = c("red", "blue")) + 
  scale_y_log10() + 
  ggtitle("casos esperados a treinta días de la primera infección")

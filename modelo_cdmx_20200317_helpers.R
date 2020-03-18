xxx <- vroom::vroom("data/Tabla_casos_positivos_resultado_InDRE_2020.03.16-Table-1-2.csv")
xxx %>% 
  mutate(date = lubridate::dmy(`Fecha de Inicio de síntomas`)) %>% 
  group_by(date) %>% 
  tally(name = "cases") %>% 
  drop_na() %>% 
  vroom::vroom_write(path = "mx.datecases.2020.03.16.tsv")


xxx %>% 
  mutate(date = lubridate::dmy(`Fecha de Inicio de síntomas`)) %>% 
  filter(Estado=="CIUDAD DE MÉXICO") %>% 
  group_by(date) %>% 
  tally(name = "cases") %>% 
  drop_na() %>% pull(cases) %>% sum
vroom::vroom_write(path = "data/cdmx.datecases.2020.03.16.tsv")

wawawa <- escenarios[,-1] %>% as.data.frame

lapply(wawawa, sum) %>% length

i = 1
parms <- as.list(scenarios[i,2:10])
init <- as.list(scenarios[i,11:28])




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

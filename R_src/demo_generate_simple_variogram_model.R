#Generates a variogram model, or adds to an existing model
vgm()
vgm("Sph")
vgm(NA, "Sph", NA, NA)
vgm(, "Sph") # "Sph" is second argument: NO nugget in this case
vgm(10, "Exp", 300)
x <- vgm(10, "Exp", 300)
vgm(10, "Nug", 0)

#sill of the model + model type + model range + nugget + smoothness parameter
vgm(10, "Exp", 300, 4.5)
vgm(10, "Mat", 300, 4.5, kappa = 0.7)
vgm( 5, "Exp", 300, add.to = vgm(5, "Exp", 60, nugget = 2.5))
vgm(10, "Exp", 300, anis = c(30, 0.5))
vgm(10, "Exp", 300, anis = c(30, 10, 0, 0.5, 0.3))
# Matern variogram model:
vgm(1, "Mat", 1, kappa=.3)
x <- vgm(0.39527463, "Sph", 953.8942, nugget = 0.06105141)
x
print(x, digits = 3);
# to see all components, do
# prints the essence of a variogram model
print.data.frame(x)
vv=vgm(model = "Tab",  covtable = 
         variogramLine(vgm(1, "Sph", 1), 1, n=1e4, min = 0, covariance = TRUE))
vgm(c("Mat", "Sph"))
vgm(, c("Mat", "Sph")) # no nugget
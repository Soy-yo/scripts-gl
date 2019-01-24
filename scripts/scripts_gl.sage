#\s
# Archivo principal que únicamente carga los demás. Cada módulo nuevo que se añada se debe cargar aquí.
# Así se pretende reducir la dependencia del directorio donde se encuentren y simplificar el uso de quienes
# tengan que copiar y pegar.
#
# Autor: Pablo Sanz Sanz
#

load("scripts/procedimientos.sage")
load("scripts/espacios.sage")
load("scripts/aplicaciones.sage")
load("scripts/recta_proyectiva.sage")
load("scripts/cuadricas.sage")
load("scripts/haces.sage")
load("scripts/afin_euclidea.sage")
load("scripts/calculadora_conicas.sage")
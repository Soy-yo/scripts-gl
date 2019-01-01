#\s
# Este archivo contiene clases para manejar cónicas y cuádricas y algunas otras funciones útiles.
# De momento no se ha implementado nada sobre cuádricas en general porque la mayoría de ejercicios son de cónicas,
# pero ciertos métodos serían exactamente iguales, sólo que en dimensión mayor. Se puede mirar el procedimiento
# y hacer a mano.
#
# Autor: Pablo Sanz Sanz
#

# Clases

#\c
# Clase que representa una cónica tanto puntual como dual dada su matriz.
#
class conica:

    #\i
    # Construye una cónica dada su matriz. Si se quiere crear mediante una ecuación es recomendable hacerlo
    # a mano (matriz simétrica con coeficientes de variables al cuadrado en la diagonal y en el resto de
    # posiciones la mitad del coeficiente según las variables). Se puede usar también ecuacion.hessian() pero
    # puede dar problemas si se tiene algún parámetro o si no aparecen las tres variables.
    #
    # Parámetros \\
    # matriz: matriz(3, 3) - matriz de la cónica
    #
    def __init__(self, matriz):
        assert matriz.is_symmetric(), "La matriz de una conica debe ser simetrica"
        assert matriz.nrows() == 3, "La matriz de una conica es 3x3"
        self._vars = vector(var('x y z'))
        self._matriz = matriz
        self._ecuacion = (self._vars * matriz * self._vars).factor()

    # Métodos accedentes

    #\m
    # Devuelve la matriz asociada a esta cónica.
    def matriz_asociada(self):
        return self._matriz

    #\m
    # Devuelve la ecuación de esta cónica.
    def ecuacion(self):
        return self._ecuacion == 0

    #\m
    # Determina si esta cónica es degenerada, esto es, si el rango de su matriz asociada no es 3.
    def es_degenerada(self):
        return self._matriz.rank() != 3

    #\m
    # Determina si esta cónica es el producto de dos rectas, esto es, si el rango de su matriz asociada es 2.
    def es_dos_rectas(self):
        return self._matriz.rank() == 2

    #\m
    # Determina si esta cónica es una recta doble, esto es, si el rango de su matriz asociada es 1.
    def es_recta_doble(self):
        return self._matriz.rank() == 1

    # Otros métodos

    #\m
    # Devuelve una nueva cónica expresada en la nueva referencia (la que sea), dada la matriz del cambio
    # de referencia.
    #
    # Implementación \\
    # Si A es la matriz asociada a esta cónica y M la matriz del cambio de referencia entre dos referencias
    # (que esta cónica desconoce), se crea una nueva cónica cuya matriz asociada sea M^t * A * M.
    #
    # Parámetros \\
    # matriz_cambio: matriz(3, 3) - matriz que representa el cambio de referencia
    #
    def cambiar_referencia(self, matriz_cambio):
        assert matriz.nrows() == 3, "La matriz del cambio de referencia de una cónica es 3x3"
        assert matriz_cambio.det() != 0, "Una matriz de cambio de referencia debe ser invertible"
        paso("La matriz de la conica en la nueva referencia se calcula:", matriz_cambio.T, self._matriz, matriz_cambio)
        return conica(matriz_cambio.T + self._matriz * matriz_cambio)

    #\m
    # Devuelve la cónica dual (de rectas tangentes) a esta cónica.
    #
    # Implementación \\
    # Se crea una nueva cónica cuya matriz asociada sea la adjunta de la de esta.
    #
    def dual(self):
        paso("La conica dual se obtiene mediante la matriz adjunta:", self._matriz.adjoint())
        return conica(self._matriz.adjoint())

    #\m
    # Devuelve la recta polar del punto P respecto de esta cónica.
    #
    # Implementación \\
    # Si A es la matriz asociada a esta cónica, la recta polar es la de ecuación x^t * A * p = 0.
    #
    # Parámetros \\
    # p: vector(3) - punto del que se quiere calcular su polar
    #
    def polar(self, p):
        assert len(p) == 3, "El punto debe ser del plano"
        paso("La ecuacion de la recta polar a p es:", matrix([self._vars]), self._matriz, matrix([p]).T, " = 0")
        v = self._matriz * p
        paso(self._vars * v, " = 0")
        # Aquí lo tenemos que hacer con dual
        _no_pasos()
        res = subespacio(v).dual()
        _no_pasos(False)
        return res

    #\m
    # Devuelve el polo de la recta r respecto de esta cónica.
    #
    # Implementación \\
    # Se calcula la recta polar del punto dual de la recta r* (el vector de sus coeficientes) respecto a la cónica dual C*.
    # Por tanto, el dual de la recta polar calculada será el punto buscado.
    #
    # Parámetros \\
    # r: subespacio - recta de la que calcular el polo
    #
    def polo(self, r):
        assert r.dim() == 1, "El polo en una conica es para rectas"
        assert r.dimension_ambiente() == 2, "La recta de la que se quiere calcular el polo debe ser recta del plano"
        _no_pasos()
        p = r.dual().representantes()[0]
        d = self.dual()
        _no_pasos(False)
        paso("Obtenemos los duales de la recta y la conica y calculamos la polar de r* (o p):")
        paso("r* = ", p, "; C*: ", d._matriz)
        polar = d.polar(p)
        _no_pasos()
        res = polar.dual().representantes()[0]
        _no_pasos(False)
        paso("Finalmente, calculamos el dual de la polar hallada")
        return res

    #\m
    # NO FUNCIONA DEL TODO: PUEDE LLAMAR A UN MÉTODO QUE AÚN NO ESTÁ IMPLEMENTADO.
    #
    # Devuelve la(s) recta(s) tangentes a la cónica que pasan por P. Será una si P pertenece a la cónica
    # y dos si no. Pueden devolverse también como cónica si así se especifica.
    #
    # Implementación \\
    # Si el punto pertenece a esta cónica devuelve su recta polar y, si no, devuelve las dos rectas que forman
    # la cónica de matriz A * p * p^t * A - (p^t * A * p) * A, donde A es la matriz asociada a esta cónica.
    # Si se especifica que se quiere oomo cónica, en cualquier caso devuelve lo segundo.
    #
    # Parámetros \\
    # p: vector(3) - punto del que obtener las tangentes \\
    # conica: booleano - determina si este método devolverá un cónica o subespacios (False (subespacios) por defecto)
    #
    def tangentes(self, p, conica = False):
        if conica:
            paso("La conica degenerada tangente a esta que pasa por el punto es la de matriz")
            paso(self._matriz, matrix([p]).T, matrix([p]), self._matriz, " - (", matrix([p]), self._matriz, \
                    matrix([p]).T, ")", self._matriz)
            return conica(self._matriz * p * p * self._matriz - (p * self._matriz * p)[0][0] * self._matriz)
        if p in self:
            paso("Como el punto pertenece a la conica, la recta tangente coincide con la polar")
            return self.polar(p)
        paso("Las rectas tangentes a esta conica que pasan por el punto son las de la conica de matriz")
        paso(self._matriz, matrix([p]).T, matrix([p]), self._matriz, " - (", matrix([p]), self._matriz, matrix([p]).T, \
                ")", self._matriz)
        c = conica(self._matriz * p * p * self._matriz - (p * self._matriz * p)[0][0] * self._matriz)
        paso("La matriz de la conica es:", c.matriz_asociada(), "; descomponemos en dos rectas")
        return c.factorizacion()

    #\m
    # Devuelve la intersección de las rectas de esta cónica degenerada. Devuelve un subespacio, si sólo
    # se quiere un punto, en caso de ser dos rectas distintas, se puede usar .representantes()[0] para
    # obtenerlo como vector. Si la cónica no es degenerada devuelve un subespacio vacío.
    #
    # Implementación \\
    # Calcula el núcleo de la matriz asociada, que será la intersección de las rectas.
    #
    def interseccion(self):
        paso("La interseccion de las rectas es el nucleo de la matriz asociada")
        paso("Usaremos el centro de una aplicacion_proyectiva por simplicidad para el codigo")
        # Aprovechamos aplicacion_proyectiva
        return aplicacion_proyectiva(self._matriz).centro()

    #\m
    # NO IMPLEMENTADO AÚN (ver pág 61).
    #
    # Devuelve una tupla conteniendo las dos rectas (como subespacios) en que se descompone esta cónica degenerada.
    # Si la cónica es una recta doble ambas rectas serán la misma.
    #
    def factorizacion(self):
        # TODO
        assert False, "No implementado"

    #\m
    # Operador in. Determina si un punto está contenido en esta cónica o no.
    #
    # Uso: P in c (P es un punto y c una cónica).
    #
    # Implementación \\
    # Comprueba p^t * A * p = 0 para A la matriz asociada de la cónica.
    #
    # Parámetros \\
    # punto: vector(3) - punto que comprobar si pertenece a la cónica
    #
    def __contains__(self, punto):
        return punto * self._matriz * punto == 0

    def __repr__(self):
        return "<Conica de ecuacion " + str(self.ecuacion()) + ">"

#\s
# Este archivo contiene clases para manejar cónicas y cuádricas y algunas otras funciones útiles.
# De momento no se ha implementado nada sobre cuádricas en general porque la mayoría de ejercicios son de cónicas,
# pero ciertos métodos serían exactamente iguales, sólo que en dimensión mayor. Se puede mirar el procedimiento
# y hacer a mano.
#
# Autor: Pablo Sanz Sanz
#

# Funciones globales

#\f
# Devuelve la cónica degenerada que contiene las dos rectas dadas. Realmente, es más útil hacerlo a mano, así que no se dará
# implementación. Es útil para otras funciones.
#
# Parámetros \\
# r1: subespacio - primera recta de la cónica degenerada \\
# r2: subespacio - segunda recta de la cónica degenerada (r1 por defecto, es decir, recta doble)
#
def conica_degenerada(r1, r2 = None):
    if r2 is None:
        r2 = r1
    assert r1.dim() == 1 and r2.dim() == 1, "Una conica degenerada esta formada por rectas"
    assert r1.dimension_ambiente() == 2 and r2.dimension_ambiente() == 2, "Las rectas deben pertenecer al plano"
    _no_pasos()
    (a1, b1, c1) = r1.dual().punto()
    (a2, b2, c2) = r2.dual().punto()
    _no_pasos(False)
    m = matrix([[2*a1*a2, a1*b2 + a2*b1, a1*c2 + a2*c1], [a1*b2 + a2*b1, 2*b1*b2, b1*c2 + b2*c1], [a1*c2 + a2*c1, b1*c2 + b2*c1, 2*c1*c2]])
    return conica(m)

#\f
# Devuelve la cónica que pasa por los 5 puntos dados.
#
# Implementación \\
# Crea el haz de cónicas que pasa por las cónicas degeneradas AB * CD y AC * BD. Después se fuerza que pase por E.
#
# Parámetros \\
# a: vector(3) - primer punto de la cónica \\
# b: vector(3) - segundo punto de la cónica \\
# c: vector(3) - tercer punto de la cónica \\
# d: vector(3) - cuarto punto de la cónica \\
# e: vector(3) - quinto punto de la cónica
#
def conica_cinco_puntos(a, b, c, d, e):
    assert len(a) == 3 and len(b) == 3 and len(c) == 3 and len(d) == 3 and len(e) == 3, "Los puntos deben pertenecer al plano"
    _no_pasos()
    ab = subespacio(a, b)
    cd = subespacio(c, d)
    abcd = conica_degenerada(ab, cd)
    ac = subespacio(a, c)
    bd = subespacio(b, d)
    acbd = conica_degenerada(ac, bd)
    _no_pasos(False)
    haz = haz_conicas(abcd, acbd)
    paso("Creamos las conicas degeneradas AB*CD y AC*BD y generamos el haz:")
    paso("(", ab.implicitas()[0], ")(", cd.implicitas()[0], ") + lambda * (", ac.implicitas()[0], ")(", bd.implicitas()[0], ")")
    paso("Por último, forzamos que pase por E =", e)
    _no_pasos()
    res = haz.forzar_punto(e)
    _no_pasos(False)
    return res

# Clases

#\c
# Clase que representa una cónica tanto puntual como dual dada su matriz. Probablemente al mostrarse su ecuación por pantalla aparezca factorizada,
# pero eso no sirve como argumento de la factorizaión.
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
    # Determina si la recta dada es tangente a esta cónica.
    #
    # Implementación \\
    # Comprueba que el punto dual de la recta pertenezca a la cónica dual.
    #
    # Parámetros \\
    # r: subespacio - recta a comprobar si es tangente o no
    #
    def es_tangente(self, r):
        assert r.dim() == 1, "El subespacio a comprobar si es tangente debe ser una recta"
        assert r.dimension_ambiente() == 2, "La recta debe pertenecer al plano"
        paso("Obtenemos la conica dual, formada por las rectas tangentes a la puntual, y el dual de la recta")
        d = self.dual()
        _no_pasos()
        p = r.dual().representantes()[0]
        _no_pasos(False)
        paso("El dual de la recta es:", p, "; comprobamos si pertenece a la conica dual")
        return p in d

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
        assert matriz_cambio.nrows() == 3, "La matriz del cambio de referencia de una cónica es 3x3"
        assert matriz_cambio.det() != 0, "Una matriz de cambio de referencia debe ser invertible"
        paso("La matriz de la conica en la nueva referencia se calcula:", matriz_cambio.T, self._matriz, matriz_cambio)
        return conica(matriz_cambio.T * self._matriz * matriz_cambio)

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
    # Devuelve una tupla conteniendo las dos rectas (como subespacios) en que se descompone esta cónica degenerada.
    # Si la cónica es una recta doble ambas rectas serán la misma. Importante ver la implementación para entender lo que se explcia en el procedimiento.
    #
    # NOTA. Aquí se utiliza el procedimiento que se describe a continuación porque es más sencillo de programar, pero hay otros métodos, como el de resolver
    # la ecuación en alguna de las variables, o hallar el punto de intersección de las rectas y luego cortar la cónica con rectas sencillas, como x=0, y=0 o z=0
    # para hallar los otros dos puntos y luego unirlos con el primero.
    #
    # Implementación \\
    # Utiliza los ejercicios 3 y 4 de la hoja 6 (han podido cambiar de numeración). Llamaremos u y v a los vectores de coeficientes de las rectas.
    # Los ejercicios nos dicen que (uv^t + vu^t)* = (uxv)(uxv)^t; uv^t - vu^t es una matriz antisimétrica con sus componentes fuera de la diagonal son, por
    # orden, (uxv)3, -(uxv)2, (uxv)1 (donde el número denota la componente del vector); y 2uv^t = (uv^t + vu^t) + (uv^t - vu^t). \\
    # Entonces, el primero de esa suma se obtiene de la propia matriz de la cónica. Para el segundo hay que obtener uxv. Eso se hace sabiendo que (uxv)(uxv)^t
    # es la matriz de la cónica dual. Calculando, la diagonal de esta matriz son los cuadrados de las componentes de uxv. Lo único que queda es calcular el signo.
    #
    def factorizacion(self):
        assert self.es_degenerada(), "La factorizacion se hace sobre una conica degenerada"
        assert self._matriz != 0, "No tiene sentido factorizar la conica de matriz 0"
        paso("Vamos a factorizar mediante un metodo visto en un ejercicio")
        m = -self._matriz.adjoint()
        paso("Primero calculamos (uxv)(uxv)^t = -A* = -", m)
        uxv = vector(self.__calcular_uxv(m))
        paso("Su daigonal principal son los cuadrados de cada componente del vector uxv => uxv = ", uxv)
        resta = matrix([[0, uxv[2], -uxv[1]], [-uxv[2], 0, uxv[0]], [uxv[1], -uxv[0], 0]])
        paso("Sabemos que uv^t-vu^t se expresa como matriz antisimetrica formada por las componentes de uxv:", resta)
        uv2 = resta + self._matriz
        paso("Y 2uv^t = ", resta, " + ", self._matriz, " = ", uv2, "; de aqui tomamos una fila y una columna y ya tenemos los coeficientes de las rectas buscadas")
        i = 1
        u = uv2.rows()[0]
        # Por si acaso
        while u == 0:
            u = uv2.rows()[i]
            i = i + 1
        i = 1
        v = uv2.columns()[0]
        # Buscamos uno linealmente independiente
        while matrix([u.list(), v.list()]).rank() == 1 and i < 3:
            v = uv2.columns()[i]
            i = i + 1
        # Era una recta doble
        if matrix([u.list(), v.list()]).rank() == 1:
            paso("Tenemos una recta doble con coeficientes: ", u)
        else:
            paso("Tenemos dos rectas de coeficientes: ", u, ", ", v)
        _no_pasos()
        res = (subespacio(u).dual(), subespacio(v).dual())
        _no_pasos(False)
        return res

    #\m
    # Calcula la imagen de un punto mediante la forma bilineal que define esta cónica. Realmente esta no es una operación de
    # cónicas, pero se incluye aquí porque puede ser útil.
    #
    # Uso: c(p) (donde c es una cónica y p un punto).
    #
    # Implementación \\
    # Calcula p^t * A * p, donde A es la matriz asociada a esta cónica.
    #
    # Parámetros \\
    # p: vector(3) - punto del que se quiere calcular su imagen
    #
    def __call__(self, p):
        return p * self._matriz * p

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
        return self(punto) == 0

    # Métodos auxiliares

    def __calcular_uxv(self, m):
        # Simplifica por si acaso
        uxv1 = sqrt(m[0][0]).simplify_full()
        uxv2 = sqrt(m[1][1]).simplify_full()
        uxv3 = sqrt(m[2][2]).simplify_full()
        # Entonces 
        if uxv1 != 0:
            # Alguno tiene signo distinto (?)
            if m[0][1] >= 0 and m[0][2] < 0:
                uxv3 = -uxv3
            elif m[0][1] < 0 and m[0][2] >= 0:
                uxv2 = -uxv2
            elif m[0][1] < 0 and m[0][2] < 0:
                uxv1 = -uxv1
        elif m[1][2] < 0:
            uxv3 = -uxv3
        return (uxv1, uxv2, uxv3)

    def __repr__(self):
        return "<Conica de ecuacion " + str(self.ecuacion()) + ">"

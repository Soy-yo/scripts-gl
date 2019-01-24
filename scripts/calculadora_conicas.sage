#\s
# Archivo que proporciona una clase que se encargará de hacer todos los cálculos de simplificación
# de una cónica, ejes, focos, etc., para un espacio euclídeo.
#
# Autor: Pablo Sanz Sanz
#

#\c
# Clase que calcula todos los elementos de una cónica calculando las referencias adecuadas, etc.
# Sus métodos deben ser llamados por orden, pues usará la información de unos para los otros.
# A pesar de los cambios de referencia que haga intermedios, los resultados serán devueltos en la
# referencia inicial que se supone canónica.
#
class calculadora_conica:

    #\i
    # Construye la calculadora dados la cónica y el espacio euclídeo. SI el espacio no tiene como
    # puntos conjugados del infinito I = (1:i:0) y J = (1:-i:0) hará un cambio de referencia para ello.
    #
    # Implementación \\
    # En caso de que I y J no sean canónicas, se cambia de referencia usando la matriz del cambio dada por
    # columnas: la parte real de cada componente de I en la primera, la parte imaginaria en la segunda y
    # tres número reales arbitrarios en la tercera que mantengan la invertibilidad de la matriz.
    #
    # Parámetros \\
    # conica: conica - cónica que se va a utilizar en la calculadora \\
    # esapcio: espacio_euclideo - espacio en el que está dado la cónica (por defecto z=0 como recta del
    # infinito y I=(1:i:0)) \\
    # col: vector(3) (real) - tercera columna de la matriz del cambio de base en caso de ser necesaria (por
    # defecto ninguna); en caso de ser necesaria pero no especificada se buscará una automáticamente (e3, e2 o e1).
    #
    def __init__(self, conica, espacio = None, col = None):
        if espacio is None:
            _no_pasos()
            espacio = espacio_euclideo()
            _no_pasos(False)
        self._conica0 = conica
        self._espacio0 = espacio
        if espacio.I() != vector([1, I, 0]):
            c1 = map(lambda x: x.real(), espacio.I())
            c2 = map(lambda x: x.imag(), espacio.I())
            if col is not None and not matrix([c1, c2, col]).is_invertible():
                col = [0, 0, 1]
            c3 = list(col) if col is not None else [0, 0, 1]
            # Alguna será invertible
            if not matrix([c1, c2, c3]).is_invertible():
                c3 = [0, 1, 0] if matrix([c1, c2, [0, 1, 0]]).is_invertible() else [1, 0, 0]
            # Matriz del cambio
            m = matrix([c1, c2, c3]).T
            paso("Cambiamos de referencia la conica y el espacio con: ", m, "^-1 = ", m^-1)
            self._cambio1 = m
            self._conica1 = conica.cambiar_referencia(self._cambio1^-1)
            self._espacio1 = espacio.cambiar_referencia(self._cambio1^-1)
            # I debería ser (1:i:0)
            paso("La matriz de la conica en esta referencia es: ", self._conica1.matriz_asociada(), "; e I = ", self._espacio1.I())
        else:
            self._cambio1 = None
            self._conica1 = conica
            self._espacio1 = espacio
        self._centro = None
        self._asintotas = None
        self._focos = None
        self._ejes = None
        self._directrices = None
        self._tipo = None
        self._cambio2 = None
        self._conica2 = None

    # Métodos accedentes

    #\m
    # Devuelve la cónica que recibió la calculadora inicialmente.
    def conica_inicial(self):
        return self._conica0

    #\m
    # Devuelve el espacio euclídeo que recibió la calculadora inicialmente,
    def espacio_inicial(self):
        return self._espacio0

    #\m
    # Devuelve la cónica que usa la calculadora para los primeros cálculos. Si I = (1:i:0) inicialmente coincidirá
    # con la cónica inicial.
    def conica_intermedia(self):
        return self._conica1

    #\m
    # Devuelve el espacio euclídeo que usa la calculadora para los cálculos. Si I = (1:i:0) inicialmente coincidirá
    # con el espacio inicial.
    def espacio_final(self):
        return self._espacio1

    #\m
    # Devuelve la cónica obtenida por el método ecuacion_reducida.
    def conica_final(self):
        return self._conica2

    # Otros métodos

    #\m
    # Calcula el centro de la cónica.
    #
    # Implementación \\
    # Como aseguramos que en la referencia actual la recta del infinito es z=0, simplemente se calculará como adj(A) * (0 0 1)^t,
    # es decir, el polo de la recta z=0 para la matriz de la cónica en la referencia actual.
    #
    def centro(self):
        adj = self._conica1.matriz_asociada().adjoint()
        paso("El centro se calcula como adj(A)", matrix([[0, 0, 1]]).T, " = ", adj, matrix([[0, 0, 1]]).T)
        c = adj * vector([0, 0, 1])
        self._centro = c / gcd(c)
        if self._centro[2] == 0:
            paso("El centro no es un punto del plano, sino que es infinito, luego la conica es una parabola")
            # Tipo parábola
            self._tipo = 2
        return self.__p(self._centro, 1)

    #\m
    # Calcula los ejes de la cónica y los devuelve como una tupla. Debe haber sido calculado el centro previamente.
    #
    # IMPORTANTE. Para parábolas no está implementado todavía.
    #
    # Implementación \\
    # Hay dos casos posibles.
    # · La cónica tiene centro. Las direcciones de los ejes son los autovectores de la submatriz de la cónica (00, 01; 10, 11).
    # Los ejes serán los que pasen por el centro y por ambas direcciones. \\
    # · La cónica no tiene centro. El propio centro infinito es la dirección del eje y la dirección perpendicular es cualquier
    # fila no nula de la submatriz anterior.
    #
    def ejes(self):
        assert self._centro is not None, "Antes de calcular los ejes se debe calcular el centro"
        return self.__ejes_p() if self._centro[2] == 0 else self.__ejes()

    #\m
    # Hace un nuevo cambio de referencia euclídeo para devolver la cónica en ecuación canónica alfa*x^2 + beta*y^2 + ganma*z^2 = 0.
    #
    # IMPORTANTE. Para parábolas no está implementado todavía.
    #
    # Implementación \\
    # Hay dos casos posibles.
    # · La cónica tiene centro. Usa las direcciones de los ejes calculadas previamente, junto con el centro como triángulo autopolar
    # de referencia y como punto unidad la suma de los anteriores, con las direcciones normalizadas (u, v) y el centro con tercera
    # coordenada 1 (c). Así, alfa = u^tAu, beta = v^tAv y ganma = c^tAc. \\
    # ·
    #
    def ecuacion_canonica(self):
        assert self._ejes is not None, "Antes de calcular la ecuacion reducida se deben calcular los ejes"
        return self.__canonica_p() if self._centro[2] == 0 else self.__canonica()

    # Métodos auxiliares

    # TODO
    def __ejes_p(self):
        assert False, "No implementado aun"

    def __ejes(self):
        subm = self._conica1.matriz_asociada()[0:2,0:2]
        paso("Obtenemos las direcciones de los ejes como los autovectores de la submatriz: ", subm)
        # Uso aplicaciones porque ya está hecho allí
        (x, y) = aplicacion_proyectiva(subm).puntos_fijos()
        u = vector(x.simplify_full().list() + [0])
        v = vector(y.simplify_full().list() + [0])
        paso("Las direcciones de los ejes son: ", u, ", ", v, "; unimos con el centro: ", self._centro)
        e1 = subespacio(self._centro, u)
        e2 = subespacio(self._centro, v)
        # Guardamos eje y dirección por si nos es útil
        self._ejes = ((u, e1), (v, e2))
        return (self.__r(e1, 1), self.__r(e2, 1))

    def __canonica(self):
        (u, v) = map(lambda x: x[0].normalized(), self._ejes)
        # La tercera coordenada es no nula
        c = self._centro / self._centro[2]
        unidad = u + v + c
        paso("Usamos como referencia euclidea R={", u, ", ", v, ", ", c, "; ", unidad, "}")
        m = matrix([u, v, c]).T.simplify_full()
        paso("Obtenemos los nuevos coeficientes de la ecuacion de la conica con: ", m)
        mc = self._conica1.cambiar_referencia(m^-1).matriz_asociada()
        paso(mc)
        (alfa, beta, ganma) = mc.diagonal()
        paso("Los ajustamos un poco si es necesario para que quede todo ordenado")
        div = gcd([alfa, beta, ganma])
        (alfa, beta, ganma) = (alfa / div, beta / div, ganma / div)
        if sign(alfa) == sign(beta) and sign(alfa) == sign(ganma):
            # Tipo elipse imaginaria
            self_tipo = -1
        else:
            if sign(alfa) == sign(beta):
                # Tipo elipse real
                self._tipo = 0
            else:
                # Tipo hipérbola
                self._tipo = 1
            if abs(alfa) > abs(beta):
                paso("Si cambiamos de orden: ", u, ", ", v, " en la referencia tendremos el eje mayor horizontal")
                (alfa, beta) = (beta, alfa)
                m = matrix([v, u, c]).T.simplify_full()
                paso("Entonces la matriz del cambio hubiera sido: ", m)
            # Hacemos que sea positiva
            if alfa < 0:
                (alfa, beta, ganma) = (-alfa, -beta, -ganma)
            (alfa, beta, ganma) = (alfa / abs(ganma), beta / abs(ganma), sign(ganma))
        self._conica2 = conica(matrix([[alfa, 0, 0], [0, beta, 0], [0, 0, ganma]]))
        self._cambio2 = m^-1
        return self._conica2

    # TODO
    def __canonica_p(self):
        assert False, "No implementado aun"

    # TODO
    def __p(self, p, n):
        paso("El punto en la referencia actual es: ", p)
        if n > 1:
            pass
        if self._cambio1 is None:
            return p
        paso("Cambiamos a la referencia canonica con la matriz del cambio: ", self._cambio1)
        return self._cambio1 * p

    # TODO
    def __r(self, r, n):
        _no_pasos()
        d = r.dual().punto()
        _no_pasos(False)
        paso("La recta en la referencia actual tiene coeficientes: ", d)
        if n > 1:
            pass
        if self._cambio1 is None:
            return r
        paso("Cambiamos a la referencia canonica con la matriz del cambio traspuesta: ", self._cambio1.T, matrix([d]).T)
        _no_pasos()
        r = subespacio(self._cambio1.T * d).dual()
        _no_pasos(False)
        return r

    def __repr__(self):
        return "<Calculadora para la conica " + str(self._conica0.ecuacion()) + " en el espacio con I = " + \
                str(self._espacio0.I()) + ">"
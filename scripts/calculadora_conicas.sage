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
            paso("Cambiamos de referencia la conica y el espacio con: ", m)
            self._cambio1 = m
            self._conica1 = conica.cambiar_referencia(self._cambio1^-1)
            self._espacio1 = espacio.cambiar_referencia(self._cambio1^-1)
            # I debería ser (1:i:0)
            paso("La matriz de la conica en esta referencia es: ", self._conica1.matriz_asociada(), "; e I = ", self._espacio1.I())
        else:
            self._cambio1 = None
            self._conica1 = conica
            self._espacio1 = espacio
        # Están en la referencia segunda (en general coincidirá la inicial) salvo que se diga lo contrario
        self._centro = None
        self._ejes = None
        self._asintotas = None
        self._focos = None # REFERENCIA FINAL
        self._vertice = None
        self._directrices = None # REFERENCIA FINAL
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
    # Devuelve la primera matriz del cambio necesaria para que I = (1:i:0). En caso de no haber sido necesario tal cambio
    # devolverá la matriz identidad.
    def matriz_cambio_inicial(self):
        return self._cambio1^-1 if self._cambio1 is not None else identity_matrix(3)

    #\m
    # Devuelve la cónica que usa la calculadora para los primeros cálculos. Si I = (1:i:0) inicialmente coincidirá
    # con la cónica inicial.
    def conica_intermedia(self):
        return self._conica1

    #\m
    # Devuelve la cónica obtenida por el método ecuacion_canonica.
    def conica_final(self):
        assert self._conica2 is not None, "Aun no se ha obtenido la ecuacion canonica"
        return self._conica2

    #\m
    # Devuelve el espacio euclídeo que usa la calculadora para los cálculos. Si I = (1:i:0) inicialmente coincidirá
    # con el espacio inicial.
    def espacio_final(self):
        return self._espacio1

    #\m
    # Devuelve la segunda matriz del cambio necesaria para que la cónica tuviera ecuación canónica.
    def matriz_cambio_final(self):
        assert self._cambio2 is not None, "Aun no se ha obtenido la ecuacion canonica"
        return self._cambio2^-1

    #\m
    # Devuelve una cadena de texto representando el tipo de la cónica que usa la calculadora. En caso de no haber sido calculado
    # aún devolverá "desconocido". En caso de haber algún error devolverá "!!!".
    def clasificacion(self):
        if self._tipo is None:
            return "desconocido"
        if self._tipo == -1:
            return "elipse imaginaria"
        if self._tipo == 0:
            return "elipse"
        if self._tipo == 1:
            return "hiperbola"
        if self._tipo == 2:
            return "parabola"
        return "!!!"

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
        return self.__p(self._centro, 0)

    #\m
    # Calcula los ejes de la cónica y los devuelve como una tupla. Debe haber sido calculado el centro previamente.
    #
    # Implementación \\
    # Hay dos casos posibles.
    # · La cónica tiene centro. Las direcciones de los ejes son los autovectores de la submatriz de la cónica (00, 01; 10, 11).
    # Los ejes serán los que pasen por el centro y por ambas direcciones. \\
    # · La cónica no tiene centro. El propio centro infinito es la dirección del eje y la dirección perpendicular es cualquier
    # fila no nula de la submatriz anterior. Se devuelve el esta otra recta como el otro eje.
    #
    def ejes(self):
        assert self._centro is not None, "Antes de calcular los ejes se debe calcular el centro"
        return self.__ejes_p() if self._centro[2] == 0 else self.__ejes()

    #\m
    # Devuelve el vértice de la cónica si lo tiene y si ha sido ya calculado. En caso de no encontrarlo dará error.
    #
    def vertice(self):
        assert self._tipo is None or self._tipo == 2, "Solo las parabolas tienen vertice"
        assert self._vertice is not None, "El vertice no ha sido calculado aun; se calcula con los ejes"
        return self._vertice

    #\m
    # Hace un nuevo cambio de referencia euclídeo para devolver la cónica en ecuación canónica alfa*x^2 + beta*y^2 + ganma*z^2 = 0.
    #
    # Implementación \\
    # Hay dos casos posibles.
    # · La cónica tiene centro. Usa las direcciones de los ejes calculadas previamente, junto con el centro como triángulo autopolar
    # de referencia y como punto unidad la suma de los anteriores, con las direcciones normalizadas (u, v) y el centro con tercera
    # coordenada 1 (c). Así, alfa = u^tAu, beta = v^tAv y ganma = c^tAc y sólo queda simplificar y reordenar. \\
    # · La cónica no tiene centro. Usa las direcciones de los ejes calculados previamente, junto con el vértice como triángulo
    # autopolar de referencia y como punto unidad la suma de los anteriores, con las direcciones normalizadas (u, v) y el vértice
    # con tercera coordenada 1 (w). Así, alfa = v^tAv y beta = w^tAu y sólo queda simplificar y reordenar.
    #
    def ecuacion_canonica(self):
        assert self._ejes is not None, "Antes de calcular la ecuacion reducida se deben calcular los ejes"
        return self.__canonica_p() if self._centro[2] == 0 else self.__canonica()

    #\m
    # Calcula las asintotas de la cónica en caso de que sea una hipérbola, partiendo de su ecuación canónica. Devuelve una tupla.
    # Si se quiere una forma más directa se recomienda usar la clase espacio_afin, que no necesita de ecuaciones reducidas.
    #
    # Implementación \\
    # Estando la hipérbola en su ecuación reducida, simplemente unimos su centro con la intersección de la cónica con la recta
    # del infinito, que se puede parametrizar como (theta : 1 : 0).
    #
    def asintotas(self):
        assert self._conica2 is not None, "Antes de calcular las asintotas se debe calcular la ecuacion canonica"
        assert self._tipo == 1, "La conica debe ser una hiperbola"
        theta = var('theta', latex_name = '\\theta')
        inf = vector([theta, 1, 0])
        ec = inf * self._conica2.matriz_asociada() * inf
        paso("Intersecamos la conica con la recta del infinito parametrizada como: ", inf)
        paso(matrix([inf]), self._conica2.matriz_asociada(), matrix([inf]).T, " = ", ec," = 0")
        sol = solve(ec, theta)
        paso(sol)
        paso("Cambiamos los puntos a las coordenadas iniciales")
        u1 = self.__p(inf.substitute(sol[0]), -1)
        u2 = self.__p(inf.substitute(sol[1]), -1)
        paso("Unimos con el centro: ", self._centro)
        a1 = subespacio(self._centro, u1)
        a2 = subespacio(self._centro, u2)
        # Guardamos asíntota y dirección por si nos es útil
        self._asintotas = ((u1, a1), (u2, a2))
        return (self.__r(a1, 0), self.__r(a2, 0))

    #\m
    # Calcula los focos de la cónica una vez está ya en su ecuación canónica. Devuelve una lista con los focos reales, dos en
    # caso de elipse e hipérbola, uno en caso de parábola.
    #
    # Implementación \\
    # Hay tres casos posibles.
    # · Elipse. Son los puntos de coordenadas (+/-sqrt(a^2 - b^2), 0), si la ecuación es (X/a)^2 + (Y/b)^2 = 1.
    # · Hipérbola. Son los puntos de coordenadas (+/-sqrt(a^2 + b^2), 0), si la ecuación es (X/a)^2 - (Y/b)^2 = 1.
    # · Parábola. Es el punto (a, 0), si la ecuación es Y^2 = 4aX.
    #
    def focos(self):
        assert self._conica2 is not None, "Antes de calcular los focos se debe calcular la ecuacion canonica"
        assert self._tipo >= 0, "No se pueden calcular los focos reales de una elipse imaginaria"
        m = self._conica2.matriz_asociada()
        # Elipse
        if self._tipo == 0:
            a2 = 1 / m[0][0]
            b2 = 1 / m[1][1]
            f1 = vector([sqrt(a2 - b2), 0, 1])
            f2 = vector([-sqrt(a2 - b2), 0, 1])
            paso("a^2 = ", a2, ", b^2 = ", b2, " => F1 = ", f1, ", F2 = ", f2)
            self._focos = [f1, f2]
        # Hipérbola
        elif self._tipo == 1:
            a2 = 1 / m[0][0]
            b2 = -1 / m[1][1]
            f1 = vector([sqrt(a2 + b2), 0, 1])
            f2 = vector([-sqrt(a2 + b2), 0, 1])
            paso("a^2 = ", a2, ", b^2 = ", b2, " => F1 = ", f1, ", F2 = ", f2)
            self._focos = [f1, f2]
        # Parábola
        else:
            a = -m[0][2] / 2
            f = vector([a, 0, 1])
            paso("a = ", a, " => F = ", f)
            self._focos = [f]
        return map(lambda f: self.__p(f, 1), self._focos)

    #\m
    # Calcula la excentricidad de la cónica, esto es, sqrt(1 - b^2/a^2) para elipses, sqrt(1 + b^2/a^2) para hipérbolas y
    # 1 para parábolas.
    #
    def excentricidad(self):
        assert self._conica2 is not None, "Antes de calcular la excentricidad se debe calcular la ecuacion canonica"
        assert self._tipo >= 0, "No se puede calcular la excentricidad de una elipse imaginaria"
        m = self._conica2.matriz_asociada()
        # Elipse
        if self._tipo == 0:
            a2 = 1 / m[0][0]
            b2 = 1 / m[1][1]
            exc = sqrt(1 - b2 / a2)
            paso("a^2 = ", a2, ", b^2 = ", b2, " => e = ", exc)
            return exc
        # Hipérbola
        elif self._tipo == 1:
            a2 = 1 / m[0][0]
            b2 = -1 / m[1][1]
            exc = sqrt(1 + b2 / a2)
            paso("a^2 = ", a2, ", b^2 = ", b2, " => e = ", exc)
            return exc
        # Parábola
        else:
            paso("La excentricidad de una parabola es 1")
            return 1

    #\m
    # Calcula las directrices de la cónica, esto es, las polares de los focos, una vez está ya en su ecuación canónica.
    # Devuelve una lista con las directrices reales, dos en caso de elipse e hipérbola, una en caso de parábola.
    #
    # Implementación \\
    # Hay tres casos posibles.
    # · Cónica con centro. Son las rectas X = a/e, X = -a/e, donde e es la excentricidad y la ecuación de la cónica es
    # (X/a)^2 + (Y/b)^2 = 1.
    # · Cónica sin centro. Es la recta X = -a, donde Y^2 = 4aX es la ecuación de la parábola.
    #
    def directrices(self):
        assert self._conica2 is not None, "Antes de calcular las directrices se debe calcular la ecuacion canonica"
        assert self._tipo >= 0, "No se pueden calcular las directrices de una elipse imaginaria"
        m = self._conica2.matriz_asociada()
        # Parábola
        if self._tipo == 2:
            a = -m[0][2] / 2
            paso("a = ", a, "; la directriz es X = ", -a)
            _no_pasos()
            self._directrices = [subespacio(vector([1, 0, a])).dual()]
            _no_pasos(False)
        # Otras
        else:
            a = sqrt(1 / m[0][0])
            exc = self.excentricidad()
            paso("a = ", a, ", e = ", exc, "; las directrices son X = +/-", a, "/", exc)
            _no_pasos()
            d1 = subespacio(vector([1, -a / exc, 0])).dual()
            d2 = subespacio(vector([1, a / exc, 0])).dual()
            _no_pasos(False)
            self._directrices = [d1, d2]
        return map(lambda d: self.__r(d, 1), self._directrices)

    # Métodos auxiliares

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
        return (self.__r(e1, 0), self.__r(e2, 0))

    def __ejes_p(self):
        subm = self._conica1.matriz_asociada()[0:2,0:2]
        paso("Obtenemos la direccion de la tangente en el vertice como cualquiera de las filas no nulas de la submatriz: ", subm)
        u = self._centro
        v = vector((subm[0].list() if subm[0] != 0 else subm[1].list()) + [0])
        paso("Las direcciones de los ejes son: ", u, ", ", v, "; nos queda obtener el vertice")
        paso("El vertice es el punto de tangencia (finito) con la conica desde: ", v, "; obtenemos la conica tangente")
        tang = self._conica1.tangentes(v).matriz_asociada()
        w = vector([tang[2][0], tang[2][1], tang[2][2] / 2])
        w = w / gcd(w)
        paso("La matriz de la conica tangente resulta: ", tang, "; tomamos su ultima fila (con la ultima coordenada entre 2) como w = ", w)
        adj = self._conica1.matriz_asociada().adjoint()
        self._vertice = adj * w
        self._vertice = self._vertice / gcd(self._vertice)
        paso("El vertice se obtiene con adj(A)", matrix([w]).T, " = ", adj, matrix([w]).T, ", proporcional a ", matrix([self._vertice]).T)
        e1 = subespacio(self._vertice, u)
        e2 = subespacio(self._vertice, v)
        # Guardamos eje y dirección por si nos es útil
        self._ejes = ((u, e1), (v, e2))
        return (self.__r(e1, 0), self.__r(e2, 0))

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
                if abs(alfa) > abs(beta):
                    paso("Si cambiamos de orden: ", u, ", ", v, " en la referencia tendremos el eje mayor horizontal")
                    (alfa, beta) = (beta, alfa)
                    m = matrix([v, u, c]).T.simplify_full()
                    paso("Entonces la matriz del cambio hubiera sido: ", m)
                # Hacemos que sea positiva
                if alfa < 0:
                    (alfa, beta, ganma) = (-alfa, -beta, -ganma)
            else:
                # Tipo hipérbola
                self._tipo = 1
                # Hay que cambiar de signo y quizá de orden
                if ganma > 0:
                    (alfa, beta, ganma) = (-alfa, -beta, -ganma)
                if alfa < 0:
                    paso("Si cambiamos de orden: ", u, ", ", v, " en la referencia tendremos la coordenada x positiva")
                    (alfa, beta) = (beta, alfa)
                    m = matrix([v, u, c]).T.simplify_full()
                    paso("Entonces la matriz del cambio hubiera sido: ", m)
            (alfa, beta, ganma) = (alfa / abs(ganma), beta / abs(ganma), sign(ganma))
        self._conica2 = conica(matrix([[alfa, 0, 0], [0, beta, 0], [0, 0, ganma]]))
        self._cambio2 = m
        return self._conica2

    def __canonica_p(self):
        (u, v) = map(lambda x: x[0].normalized(), self._ejes)
        # La tercera coordenada es no nula
        w = self._vertice / self._vertice[2]
        unidad = u + v + w
        paso("Usamos como referencia euclidea R={", u, ", ", v, ", ", w, "; ", unidad, "}")
        m = matrix([u, v, w]).T.simplify_full()
        paso("Obtenemos los nuevos coeficientes de la ecuacion de la conica con: ", m)
        mc = self._conica1.cambiar_referencia(m^-1).matriz_asociada()
        paso(mc)
        (alfa, beta) = (mc[1][1], mc[0][2])
        p = beta / alfa
        paso("Los ajustamos un poco si es necesario para que quede todo ordenado")
        if p > 0:
            paso("Cambiamos la orientacion de: ", u, " en la referencia para tener la parabola hacia la derecha")
            m = matrix([-u, v, w]).T.simplify_full()
            p = -p
        self._conica2 = conica(matrix([[0, 0, p], [0, 1, 0], [p, 0, 0]]))
        self._cambio2 = m
        return self._conica2

    # n = 0 => cambio a inicial
    # n = 1 => cambio a intermedia y a inicial
    # n = -1 => cambio a intermedia
    def __p(self, p, n):
        paso("El punto en la referencia actual es: ", p)
        if abs(n) > 0 and self._cambio2 is not None:
            paso("Cambiamos a la referencia con la que estabamos trabajando con la matriz del cambio: ", self._cambio2)
            p = self._cambio2 * p
            paso(p)
        if self._cambio1 is None or n < 0:
            return p
        paso("Cambiamos a la referencia canonica con la matriz del cambio: ", self._cambio1)
        return self._cambio1 * p


    # n = 0 => cambio a inicial
    # n = 1 => cambio a intermedia y a inicial
    # n = -1 => cambio a intermedia
    def __r(self, r, n):
        _no_pasos()
        d = r.dual().punto()
        _no_pasos(False)
        paso("La recta en la referencia actual tiene coeficientes: ", d)
        if abs(n) > 0 and self._cambio2 is not None:
            paso("Cambiamos a la referencia con la que estabamos trabajando con la matriz del cambio menos traspuesta: ", self._cambio2.T^-1)
            d = self._cambio2.T^-1 * d
            paso(d)
        if self._cambio1 is None or n < 0:
            _no_pasos()
            r = subespacio(d).dual()
            _no_pasos(False)
            return r
        paso("Cambiamos a la referencia canonica con la matriz del cambio menos traspuesta: ", self._cambio1.T^-1)
        _no_pasos()
        r = subespacio(self._cambio1.T^-1 * d).dual()
        _no_pasos(False)
        return r

    def __repr__(self):
        return "<Calculadora para la conica " + str(self._conica0.ecuacion()) + " en el espacio con I = " + \
                str(self._espacio0.I()) + ">"
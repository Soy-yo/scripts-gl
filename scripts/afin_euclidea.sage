#\s
# Archivo con clases representantes de espacios tanto afines como euclídeos con diversas funciones relacionadas con
# aplicaciones, cónicas, etc.
#
# Autor: Pablo Sanz Sanz
#

# Clases

#\c
# Clase que representa un espacio afín de una dimensión arbitraria. Sirve para diversos cálculos relacionados con
# subespacios, cónicas, etc.
#
# NOTA. Se puede indicar que por defecto los resultados sean devueltos en la referencia canónica o no. En cualquier caso los
# puntos que se pasan como parámetro deben estar en la referencia de este espacio.
# Si, por ejemplo, queremos calcular el punto medio de los puntos (-1 : 0 : 1) y (1 : 0 : 1) (dados en la referencia canónica),
# pero la referencia de este espacio no es la canónica, hay que cambiar primero esos puntos a la referencia de este espacio a
# mano con la inversa de la matriz que devuelva el método matriz_cambio (M^-1 * X) y después llamar al método con estos
# dos nuevos puntos. Finalmente, en función del valor que se haya indicado a cambio_defecto, se cambiará la referencia del
# resultado no.
#
# Todo esto se aplica en general, salvo que se indique lo contrario. En el caso de métodos que devuelvan subespacios (p.e:
# asintotas_conica) lo normal será que no se cambie la referencia automáticamente por defecto. Importante tenerlo en cuenta.
#
class espacio_afin:

    #\i
    # Construye un espacio afín, dados el hiperplano del infinito y la referencia que se está utilizando.
    # Nótese que el hiperplano se supone que viene expresado en la referencia dada.
    #
    # Parámetros \\
    # hiperplano_infinito: subespacio(dim=n-1) - hiperplano que se está considerando como infinito (xn = 0 por defecto) \\
    # dimension: natural - dimensión de este espacio si no se han especificado los dos anteriores parámetros; por ejemplo,
    # espacio_afin(dimension = 2) crea un espacio afín cuya recta del infinto es z = 0 y la referencia es la canónica
    # (2 por defecto)
    #
    def __init__(self, hiperplano_infinito = None, dimension = 2):
        self._dim = hiperplano_infinito.dimension_ambiente() if hiperplano_infinito is not None else dimension
        _no_pasos()
        # xn = 0 por defecto
        if hiperplano_infinito is None:
            hiperplano_infinito = subespacio(vector(([0] * self._dim) + [1])).dual()
        _no_pasos(False)
        assert hiperplano_infinito.dimension_ambiente() == self._dim, "La dimension ambiente del hiperplano debe coincidir con la de este espacio"
        assert hiperplano_infinito.dim() == self._dim - 1, "El subespacio debe ser un hiperplano"
        self._infinito = hiperplano_infinito

    # Métodos accedentes

    #\m
    # Devuelve la dimensión de este espacio afín.
    def dim(self):
        return self._dim

    #\m
    # Devuelve el hiperplano del infinito de este espacio afín en la referencia en que está expresado el espacio.
    def hiperplano_infinito(self):
        return self._infinito

    # Otros métodos

    #\m
    # Determina si la referencia indicada es válida para es espacio, esto es, todos los vectores son de longitud adecuada,
    # son suficientes y todos menos el último son proyectivamente independientes.
    #
    # Parámetros \\
    # referencia: lista(vector(n+1))(long=n+2) - referencia a comprobar
    #
    def es_referencia_valida(self, referencia):
        if len(referencia) != self._dim + 2:
            return False
        if any(map(lambda p: len(p) != self._dim + 1, referencia)):
            return False
        if matrix(referencia).rank() != len(referencia) - 1:
            return False
        return True

    #\m
    # Determina si la referencia indicada es una referencia afín, esto es, cada punto de los n primeros pertenece al hiperplano
    # del infinito.
    #
    # Parámetros \\
    # referencia: lista(vector(n+1))(long=n+2) - referencia a comprobar
    #
    def es_referencia_afin(self, referencia):
        assert self.es_referencia_valida(referencia), "La referencia dada debe ser valida"
        for p in referencia[0 : -2]:
            if p not in self._infinito:
                return False
        return True

    #\m
    # Devuelve un nuevo espacio afín (no modifica este) que se supone expresado en una referencia R' dada la matriz del cambio
    # de la referencia actual del espacio a R'.
    #
    # Implementación \\
    # Cambia la referencia del hiperplano del infinito usando M^-t * u, donde u es el vector asociado a las coordenadas de H
    # (el vector dual al subespacio).
    #
    # Parámetros \\
    # matriz_cambio: matriz(n+1,n+1) - matriz del cambio de referencia de R a R'
    #
    def cambiar_referencia(self, matriz_cambio):
        assert matriz_cambio.is_square() and matriz_cambio.nrows() == self._dim + 1, "La matriz debe representar un cambio en este espacio"
        assert matriz_cambio.det() != 0, "Una matriz de cambio de referencia debe ser invertible"
        paso("Cambiamos el hiperplano del infinito de referencia, cambiando de referencia el punto dual que lo representa con " + \
                "la matriz del cambio inversa y traspuesta M^-t = ", matriz_cambio.T^-1)
        _no_pasos()
        dual = self._infinito.dual().punto()
        coord = matriz_cambio.T^-1 * dual
        h = subespacio(coord).dual()
        _no_pasos(False)
        paso("Sus coordenadas: ", dual, " se convierten en: ", coord, ": ", h.implicitas()[0])
        return espacio_afin(h)

    #\m
    # Determina los puntos del infinito del subespacio dado, esto es la intersección de este con el hiperplano del infinito.
    # Devuelve un objeto del tipo subespacio. Si el subespacio dado es una recta se puede tomar el único punto de intersección
    # con el método punto() de la clase subespacio.
    #
    # Parámetros \\
    # s: subespacio(dim_ambiente=n) - subespacio del que se quiere conocer los puntos del infinito
    #
    def puntos_infinitos(self, s):
        assert s.dimension_ambiente() == self._dim, "El subespacio debe ser de este espacio"
        paso("Intersecamos con el hiperplano del infinito: ", self._infinito.implicitas()[0])
        return s.interseccion(self._infinito)

    #\m
    # Determina si dos subespacios son paralelos, esto es, si la intersección con el hiperplano del infinito de uno esá contenida
    # (o es igual) en la del otro.
    #
    # Parámetros \\
    # r: subespacio(dim_ambiente=n) - uno de los subespacios que se quiere comprobar si son paralelos \\
    # s: subespacio(dim_ambiente=n) - el otro
    #
    def paralelos(self, r, s):
        # Nos quedamos en r con el de menor dimensión
        if r.dim() > s.dim():
            temp = r
            r = s
            s = temp
        inf1 = self.puntos_infinitos(r)
        inf2 = self.puntos_infinitos(s)
        # ¿Están todos los puntos representantes del pequeño en el grande?
        return all(map(lambda p: p in inf2, inf1.representantes()))

    #\m
    # Calcula la razón simple (A, B, C) de los puntos alineados dados.
    #
    # Implementación \\
    # Calcula la razón doble {Pinf, A; B, C}, donde Pinf es el punto de intersección entre la recta que forman los puntos y el
    # hiperplano del infinito (ver razon_doble (recta_proyectiva.sage)).
    #
    # Parámetros \\
    # a: vector(n+1) - primer punto de la razón simple \\
    # b: vector(n+1) - segundo punto de la razón simple \\
    # c: vector(n+1) - tercer punto de la razón simple
    #
    def razon_simple(self, a, b, c):
        assert len(a) == self._dim + 1 and len(b) == self._dim + 1 and len(c) == self._dim + 1, "Los puntos deben pertenecer a este espacio"
        paso("Obtenemos la recta por la que pasan los puntos:")
        r = subespacio(a, b, c)
        assert r.dim() == 1, "Los puntos deben estar alineados"
        pinf = self.puntos_infinitos(r).punto()
        paso("El punto del infinito es: ", pinf, "; calculamos la razon doble:")
        return razon_doble_puntos(pinf, a, b, c)

    #\m
    # Determina si los puntos dados son simétricos respecto del centro dado. Para ello deben yacer todos en la misma recta.
    #
    # Implementación \\
    # Determina si {P1, P2 ; C, Pinf} es una cuaterna armónica, donde Pinf es el punto de intersección entre la recta que forman
    # los puntos y el hiperplano del infinito (ver es_cuaterna_armonica (recta_proyectiva.sage)).
    #
    def son_simetricos(self, p1, p2, centro):
        assert len(p1) == self._dim + 1 and len(p2) == self._dim + 1 and len(centro) == self._dim + 1, "Los puntos deben pertenecer a este espacio"
        paso("Obtenemos la recta por la que pasan los puntos:")
        r = subespacio(p1, p2, centro)
        assert r.dim() == 1, "Los puntos deben estar alineados"
        pinf = self.puntos_infinitos(r).punto()
        paso("El punto del infinito es: ", pinf, "; determinamos si son cuaterna armonica:")
        return bool(es_cuaterna_armonica(p1, p2, centro, pinf))

    #\m
    # Calcula el punto simétrico a un punto respecto de un punto "centro" de la simetría. El punto calculado estará en la recta que
    # forman los puntos dados.
    #
    # Implementación \\
    # Monta la recta (ver recta_proyectiva (recta_proyectiva.sage)) que pasa por los puntos, tomando como referencia {Pinf, C; P}, donde
    # Pinf es el punto del infinito de la recta PC y C el centro dado. Entonces, el simétrico de P será aquel cuya coordenada sea la opuesta
    # a la de P, esto es, el punto de coordenada -1.
    #
    # Parámetros \\
    # punto: vector(n+1) - punto del que se quiere calcular su simétrico \\
    # centro: vector(n+1) - punto centro de la simetría
    #
    def simetrico(self, punto, centro):
        assert len(punto) == self._dim + 1 and len(centro) == self._dim + 1, "Los puntos deben pertenecer a este espacio"
        paso("Obtenemos la recta por la que pasan los puntos:")
        s = subespacio(punto, centro)
        pinf = self.puntos_infinitos(s).punto()
        paso("El punto del infinito es: ", pinf, "; montamos la recta (realmente es la misma, pero parametrizada)")
        r = recta_proyectiva(pinf, centro, punto)
        paso("El simetrico de: ", punto, " sera el punto de la recta con coordenada -1")
        return r[-1]

    #\m
    # Calcula el punto medio entre p1 y p2. Es lo contrario al método simetrico.
    #
    # Implementación \\
    # Calcula el conjugado armónico de Pinf respecto de P1 y P2, donde Pinf es el punto del infinito de la recta P1P2 (ver conjugado_armonico
    # (recta_proyectiva.sage)).
    #
    # Parámetros \\
    # p1: vector(n+1) - uno de los puntos de los que se quiere calcular el punto medio \\
    # p2: vector(n+1) - el otro punto
    #
    def punto_medio(self, p1, p2):
        assert len(p1) == self._dim + 1 and len(p2) == self._dim + 1, "Los puntos deben pertenecer a este espacio"
        paso("Obtenemos la recta por la que pasan los puntos:")
        r = subespacio(p1, p2)
        pinf = self.puntos_infinitos(r).punto()
        paso("El punto del infinito es: ", pinf, "; calculamos su conjugado armonico respecto de los otros dos:")
        return conjugado_armonico(p1, p2, pinf)

    #
    # Plano afín
    #

    #\m
    # Devuelve los dos puntos de intersección de la recta del infinito con la cónica dada.
    #
    # Parámetros \\
    # c: conica/parametrizacion_conica - cónica a intersecar
    #
    def puntos_infinitos_conica(self, c):
        assert self._dim == 2, "Las funciones sobre conicas son para el plano afin (dimension 2)"
        r = self._infinito
        if c._tipo() == 0:
            _no_pasos()
            rep = self._infinito.representantes()
            r = recta_proyectiva(rep[0], rep[1])
            _no_pasos(False)
            theta = var('theta_var', latex_name = '\\theta')
            paso("La recta del infinito se puede expresar como: ", r[theta])
        return c.interseccion_recta(r)

    #\m
    # Devuelve una cadena de texto que indica el tipo de la cónica dada por parámetro. No se diferencian casos
    # reales de imaginarios (p.e: elipse real o imaginaria).
    #
    # Implementación \\
    # Calcula la intersección de la cónica con la recta del infinito y se pueden dar los siguientes casos: \\
    # · Intersección en dos puntos reales: hipérbola. \\
    # · Intersección en dos puntos complejos: elipse \\
    # · Intersección en un punto doble: parábola
    #
    # Parámetros \\
    # c: conica/parametrizacion_conica - cónica a clasificar
    #
    def clasificacion_conica(self, c):
        t = self.__tipo_conica(c)
        return "hiperbola" if t == 0 else "elipse" if t == 1 else "parabola"

    #\m
    # Devuelve el centro de la cónica dada, esto es, el polo de la recta del infinito. Si la cónica es una parábola no tendrá centro y
    # devolverá un vector vacío.
    #
    # Parámetros \\
    # c: conica - cónica de la que se quiere saber el centro
    #
    def centro_conica(self, c):
        assert self._dim == 2, "Las funciones sobre conicas son para el plano afin (dimension 2)"
        paso("El centro es el polo de la recta del infinito")
        centro = c.polo(self._infinito)
        return centro if centro not in self._infinito else vector([])

    #\m
    # Devuelve una tupla conteniendo las asíntotas de la cónica dada. Nótese que en caso de ser una elipse estas serán rectas imaginarias
    # y en caso de ser una parábola no tendrá (pues no tiene centro) y se devolverán dos subespacios vacíos.
    #
    # Implementación \\
    # Las asíntotas son tangentes a la cónica en la recta del infinito. Se calcula la intersección con ella y el resultado serán las
    # tangentes (las polares) de cada punto.
    #
    # Parámetros \\
    # c: conica - cónica de la que obtener las asíntotas
    #
    def asintotas_conica(self, c):
        paso("Intersecamos la conica con la recta del infinito")
        (a, b) = self.puntos_infinitos_conica(c)
        if matrix([a, b]).rank() == 1:
            paso("La interseccion es un punto doble: ", a, "; no hay asintotas")
            # Subespacio vacío
            return (subespacio(), subespacio())
        paso("Calculamos las tangentes (las polares) de los puntos obtenidos: ", (a, b))
        r1 = c.polar(a)
        r2 = c.polar(b)
        return (r1, r2)

    #\m
    # Devuelve una lista con los vectores de la referencia afín en que la cónica dada tendrá ecuación reducid x^2 +/- y^2 +/- z^2 = 0.
    #
    # NOTA. Cuidado de no usar asíntotas como diámetro en el parámetro, pues son tangentes a la cónica.
    #
    # Implementación \\
    # Si la cónica tiene centro, se toma como origen de coordenadas, esto es, X2 = C. Los otros dos puntos, X0 y X1, se toman como los
    # puntos del infinito de dos diámetros conjugados cualesquiera (se debe indicar uno de ellos como parámetro). \\
    # En el caso de la parábola se tomará como origen un punto de la cónica del que se podrá indicar su coordenada x (z será 1 por defecto). Se
    # tomarán, además, X0 el punto de tangencia con la recta del infinito y el polo de X0X2 como X1. \\
    # El punto unidad será escogido en ambos casos por el método referencia_autopolar de la cónica (conica (cuadricas.sage)).
    #
    # Parámetros \\
    # c: conica - cónica de la que obtener la referencia en la que tiene ecuación reducida \\
    # diam: subespacio(dim=1,dim_amb=2) - diámetro que se quiera usar para el primer caso (por defecto ninguno) \\
    # x2: vector(3) - tercer punto que se quiera tomar en la referencia en el segundo caso (por defecto ninguno) \\
    # e: vector(3) - punto unidad que se quiera tomar en la referencia en el segundo caso (por defecto ninguno) \\
    # real: booleano - determina si la cónica debe ser necesariamente real o no (True por defecto)
    #
    def referencia_ecuacion_reducida(self, c, diam = None, x2 = None, e = None, real = True):
        paso("Calculamos el centro de la conica si lo tiene")
        centro = self.centro_conica(c)
        # Parábola
        if len(centro) == 0:
            assert x2 is not None and e is not None, "Se deben indicar los puntos X2 y unidad para la parabola"
            _no_pasos()
            # Tomamos el primero de los puntos del infinito (son el mismo)
            x_0 = self.puntos_infinitos_conica(c)[0]
            _no_pasos(False)
            paso("No hay centro, usamos el punto de tangencia que hubiera sido el centro como X0 = ", x_0)
            x_2 = x2
            paso("Usamos como X2 y unidad los puntos indicados: X2 = ", x_2, "; E = ", e)
            _no_pasos()
            r = subespacio(x_0, x_2)
            _no_pasos(False)
            paso("X1 es el polo de la recta X0X2: ", r.implicitas()[0])
            x_1 = c.polo(r)
            paso("X1 = ", x_1)
            return [x_0, x_1, x_2, e]
        # Elipse/hipérbola
        else:
            x_2 = centro
            paso("Tomamos el centro obtenido como origen de coordenadas: X2 = ", x_2)
            _no_pasos()
            b = d is not None and self.es_diametro(c, diam)
            _no_pasos(False)
            assert b, "Se debe indicar un diametro de la conica con el que continuar"
            paso("Obtenemos el diametro donjugado del dado:")
            conj = self.diametro_conjugado(c, diam)
            _no_pasos()
            x_0 = self.puntos_infinitos(diam).punto()
            x_1 = self.puntos_infinitos(conj).punto()
            _no_pasos(False)
            paso("El conjugado es: ", conj.implicitas()[0], "; y los puntos del infinito de cada diametro: X0 = ", x_0, ", X1 = ", x_1)
            paso("Completamos la referencia con el punto unidad adecuado:")
            return [x_0, x_1, x_2, c.referencia_autopolar(x_0, x_1, x_2, real)]

    #\m
    # Determina si una recta d dada como subespacio es un diámetro de la cónica indicada, es decir, si pasa por su centro.
    #
    # Parámetros \\
    # c: conica - cónica de la que se quiere saber si d es diámetro \\
    # d: subespacio(dim=1,dim_amb=2) - recta candidata a ser diámetro
    #
    def es_diametro(self, c, d):
        assert self._dim == 2, "Las funciones sobre conicas son para el plano afin (dimension 2)"
        assert d.dim() == 1, "El diametro debe ser una recta"
        assert d.dimension_ambiente() == 2, "La recta debe ser del plano"
        paso("Obtenemos el centro")
        centro = self.centro_conica(c)
        # No tiene centro
        if len(centro) == 0:
            paso("La conica no tiene centro")
            return False
        paso("El centro es: ", centro, "; comprobamos si pertenece al diametro")
        return centro in d

    #\m
    # Calcula el diámetro conjugado del dado respecto de la conica indicada.
    #
    # Implementación \\
    # Obtiene la recta polar del punto del infinito del diámetro respecto de la cónica.
    #
    # Parámetros \\
    # c: conica - cónica sobre la que se quiere calcular el diámetro conjugado \\
    # d: subespacio(dim=1,dim_amb=2) - diámetro del que se quiere calcular su conjugado
    #
    def diametro_conjugado(self, c, d):
        _no_pasos()
        b = self.es_diametro(c, d)
        _no_pasos(False)
        assert b, "La recta dada debe ser diametro de la conica"
        pinf = self.puntos_infinitos(d).punto()
        paso("El punto del infinito es: ", pinf, "; calculamos su polar")
        return c.polar(pinf)

    #\m
    # Determina si los diámetros d1 y d2 son conjugados respecto a la cónica c. Esto ocurre cuando sus puntos del infinito
    # son conjugados respecto a c, o lo que es lo mismo, que la polar del punto del infinito de uno sea el otro diámetro.
    #
    # Parámetros \\
    # c: conica - cónica de la que se quiere saber si los diámetros son conjugados \\
    # d1: subespacio(dim=1,dim_amb=2) - primer diámetro \\
    # d2: subespacio(dim=1,dim_amb=2) - segundo diámetro
    #
    def son_diametros_conjugados(self, c, d1, d2):
        _no_pasos()
        b = self.es_diametro(c, d1) and self.es_diametro(c, d2)
        _no_pasos(False)
        assert b, "Las rectas dadas deben ser diametros de la conica"
        paso("Calculamos el diametro conjugado del primero y vemos si coincide con el segundo")
        conj = self.diametro_conjugado(c, d1)
        paso("La polar tiene ecuacion: ", conj.implicitas()[0], "; el otro diametro: ", d2.implicitas()[0], "; comprobamos que coincidan")
        return conj == d2

    # Métodos auxiliares

    def __tipo_conica(self, c):
        paso("Intersecamos la conica con la recta del infinito")
        (a, b) = self.puntos_infinitos_conica(c)
        paso("La interseccion es: ", a, ", ", b)
        # Los puntos son distintos (con que uno sea real los dos deberían serlo pero se comprueba de todos modos)
        if matrix([a, b]).rank() == 2:
            return 0 if self.__es_real(a) and self.__es_real(b) else 1
        return 2

    def __es_real(self, p):
        return p == conjugate(p)

    def __repr__(self):
        return "<Espacio afin de dimension " + str(self._dim) + " con hiperplano del infinito de ecuacion " + str(self._infinito.implicitas()[0]) + ">"
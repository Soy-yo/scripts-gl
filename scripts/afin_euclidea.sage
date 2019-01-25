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
class espacio_afin:

    #\i
    # Construye un espacio afín, dado el hiperplano del infinito.
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
    # de la referencia canónica a R'.
    #
    # Implementación \\
    # Cambia la referencia del hiperplano del infinito usando M^-t * u, donde u es el vector asociado a las coordenadas de H
    # (el vector dual al subespacio).
    #
    # Parámetros \\
    # matriz_cambio: matriz(n+1,n+1) - matriz del cambio de referencia a R'
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
    # Devuelve la dirección de la recta indicada, esto es, su corte con la recta del infinito.
    #
    # Parámetros \\
    # r: subespacio(dim=1,dim_amb=n) - recta de la que calcular su dirección
    #
    def direccion(self, r):
        assert r.dim() == 1, "El subespacio debe ser una recta"
        return self.puntos_infinitos(r).punto()

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
    # Calcula la recta paralela la recta r que pasa por el punto P.
    #
    # Implementación \\
    # Calcula el punto del infinito de la recta y lo une con P.
    #
    # Parámetros \\
    # r: subespacio(dim=1,dim_amb=n) - recta de la que obtener la paralela \\
    # p: vector(m+1) - punto por el que debe pasar la paralela
    #
    def paralela(self, r, p):
        assert len(p) == self._dim + 1, "El punto debe pertenecer a este espacio"
        paso("Calculamos la direccion de la recta dada, que es la interseccion con el hiperplano del infinito")
        _no_pasos()
        u = self.direccion(r)
        _no_pasos(False)
        paso(u, "; unimos con el punto dado: ", p)
        return subespacio(p, u)

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
    # Devuelve el punto que se obtiene al desplazar el punto P dado en la dirección del vector AB y distancia |AB|.
    #
    # Implementación \\
    # Calcula la recta paralela desde P al vector AB y halla la intersección con la paralela a PA desde B.
    #
    # Parámetros \\
    # p: vector(n+1) - punto a desplazar \\
    # a: vector(n+1) - punto de origen del vector \\
    # b: vector(n+1) - punto de destino del vector
    #
    def desplazar(self, p, a, b):
        assert len(p) == self._dim + 1 and len(a) == self._dim + 1 and len(b) == self._dim + 1, "Los puntos deben pertenecer a este espacio"
        paso("Calculamos la recta que une los puntos del vector: ", a, ", ", b , ", para luego calcular la paralela a esta que pasa por el punto dado:")
        _no_pasos()
        ab = subespacio(a, b)
        _no_pasos(False)
        paso("La recta es: ", ab.implicitas()[0])
        r = self.paralela(ab, p)
        paso("Calculamos la recta que une el punto: ", p, " con el primer punto del vector: ", a, ", para luego calcular su paralela por: ", b)
        _no_pasos()
        pa = subespacio(p, a)
        _no_pasos(False)
        paso("La recta es: ", pa.implicitas()[0])
        s = self.paralela(pa, b)
        paso("El resultado es la interseccion")
        return r.interseccion(s).punto()

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
            return 0 if es_real(a) and es_real(b) else 1
        return 2

    def __repr__(self):
        return "<Espacio afin de dimension " + str(self._dim) + " con hiperplano del infinito de ecuacion " + str(self._infinito.implicitas()[0]) + ">"

#\c
# Clase que representa un espacio euclídeo de dimensión 2. Sirve para diversos cálculos relacionados con subespacios, cónicas, etc.
# No utiliza la cuádrica del absoluto, pero sí los puntos conjugados del infinito.
#
# IMPORTANTE. Cuidado de no llamar a ninguna variable I por ningún lugar o probablemente esto deje de funcionar. \\
# IMPORTTANTE 2. A pesar de que esta clase contiene métodos para cálculo de elementos de cónicas, utilizan un método "directo", desde la definición,
# por lo tanto, los resultados obtenidos pueden no ser correctos. En caso de querer calcular focos, ejes, etc., se recomienda mucho más usar
# la clase calculadora_conica (calculadora_conica.sage), que realiza los cálculos simplificando la cónica primero, con lo cual sí deberían estar
# bien en cualquier caso.
#
# NOTA. Esta clase está basada en espacio_afin, pero no me he molestado en aprender herencia en Python, así que, si se quieren usar métodos de
# espacio_afin, se deberá acceder al que esta clase proporciona mediante el método espacio_afin().
#
class espacio_euclideo:

    #\i
    # Construye un espacio euclídeo dado uno de sus puntos conjugados del infinito.
    #
    # Parámetros \\
    # i: vector(3)(imaginario) - uno de los puntos cíclicos conjugados del infinito, del que se obtiene la recta del infinito ((1 : i : 0) por defecto)
    #
    def __init__(self, i = vector([1, I, 0])):
        assert len(i) == 3, "El punto conjugado del infinito debe pertenecer al plano"
        self._i = i
        # Por vaguería
        self._j = conjugate(i)
        paso("Obtenemos la recta del infinito a partir de los puntos conjugados del infinito I y J")
        infinito = subespacio(self._i, self._j)
        self._afin = espacio_afin(infinito)

    # Métodos accedentes

    #\m
    # Devuelve el espacio afín asociado a este espacio euclídeo.
    def espacio_afin(self):
        return self._afin

    #\m
    # Devuelve la recta del infinito de este espacio euclídeo.
    def recta_infinito(self):
        return self._afin.hiperplano_infinito()

    #\m
    # Devuelve una tupla conteniendo los puntos conjugados del infinito.
    def conjugados_infinito(self):
        return (self._i, self._j)

    #\m
    # Devuelve el primer punto conjugado del infinito.
    def I(self):
        return self._i

    #\m
    # Devuelve el segundo punto conjugado del infinito.
    def J(self):
        return self._j

    # Otros métodos

    #\m
    # Devuelve un nuevo espacio euclídeo (no modifica este) que se supone expresado en una referencia R' dada la matriz del cambio
    # de la referencia canónica a R'.
    #
    # Implementación \\
    # Cambia la referencia del espacio afín asociado a este y las coordenadas de los puntos conjugados del infinito.
    #
    # Parámetros \\
    # matriz_cambio: matriz(3,3) - matriz del cambio de referencia a R'
    #
    def cambiar_referencia(self, matriz_cambio):
        i2 = matriz_cambio * self._i
        paso("Calculamos las nuevas coordenadas de I: ", matriz_cambio, matrix([self._i]).T, " = ", i2)
        return espacio_euclideo(i2)

    #\m
    # Determina si las direcciones son ortogonales.
    #
    # Para determinar si lo son dos rectas, primero calcular su dirección.
    #
    # Implementación \\
    # Determina si los puntos forman una cuaterna armónica con los puntos conjugados del infinito.
    #
    # u: vector(3) - uno de los puntos del infinito de los que se quiere saber si son perpendiculares \\
    # v: vector(3) - el otro
    #
    def perpendiculares(self, u, v):
        assert u in self.recta_infinito() and v in self.recta_infinito(), "Los puntos deben ser puntos del infinito"
        paso("Dos direcciones son perpendiculares si forman una cuaterna armonica con los conjugados del infinito: ", self._i, ", ", self._j)
        return bool(es_cuaterna_armonica_puntos(u, v, self._i, self._j))

    #\m
    # Determina la dirección perpendicular a la dada.
    #
    # Implementación \\
    # Calcula el cuarto armónico de u respecto de los puntos conjugados del infinito.
    #
    # Parámetros \\
    # u: vector(3) - el punto del infinito del que se quiere saber su perpendicular
    #
    def direccion_perpendicular(self, u):
        assert u in self.recta_infinito(), "El punto debe ser del infinito"
        paso("La direccion perpendicular sera el cuarto armonico con los conjugados del infinito: ", self._i, ", ", self._j)
        return conjugado_armonico_puntos(self._i, self._j, u)

    #\m
    # Determina la recta perpendicular a la dada que pasa por el punto P dado.
    #
    # Implementación \\
    # Calcula la dirección perpendicular y une con P.
    #
    # Parámetros \\
    # u: vector(3) - el punto del infinito del que se quiere saber su perpendicular \\
    # p: vector(3) - punto por el que debe pasar la perpendicular calculada
    #
    def perpendicular(self, u, p):
        v = self.direccion_perpendicular(u)
        return subespacio(p, v)

    #\m
    # Calcula el ángulo entre dos direcciones.
    #
    # Para calcularlo para rectas, calcular su dirección primero.
    #
    # Implementación \\
    # Calcula 1/2i * Log{u, v; I, J}
    #
    # u: vector(3) - uno de los puntos del infinito de los que se quiere saber si son perpendiculares \\
    # v: vector(3) - el otro
    #
    def angulo(self, u, v):
        paso("La formula para el angulo es 1/2i * Log{", u, ", ", v, "; ", self._i, ", ", self._j, "}")
        r = razon_doble_puntos(u, v, self._i, self._j)
        paso("La razon doble es: ", r)
        return 1/(2*I) * log(r)

    #\m
    # Determina si la cónica C dada es una circunferencia, esto es, si contiene los puntos conjugados del infinito.
    #
    # Parámetros \\
    # c: conica/parametrizacion_conica - conica a comprobar si es una circunferencia
    #
    def es_circunferencia(self, c):
        return self._i in c and self._j in c

    #\m
    # Devuelve una lista contieniendo los focos de la cónica dada.
    #
    # NOTA. Puede ser un poco lento. Paciencia. \\
    # IMPORTANTE. No es fiable 100% cuando depende de la factorización de cónicas degeneradas, pues no es un método fiable 100%.
    # Si se obtienen resultados muy raros o errores es recomendable ejecutar varias veces. Si eso no lo arregla, tocará hacerlo a mano
    # (seguir la implementación que se da más abajo, pero hallando las rectas tangentes a mano).
    #
    # Implementación \\
    # Se obtienen las cónicas degeneradas tangentes a c desde I y J y se obtienen los puntos de corte F, F', G y G', donde G, G'
    # son los focos imaginarios. En caso de ser una circunferencia, únicamente devuelve su centro.
    #
    # Parámetros \\
    # c: conica - cónica de la que calcular los focos
    #
    def focos(self, c):
        if S.es_circunferencia(c):
            paso("I = ", self._i, ", J = ", self._j, " pertenecen a la conica, luego es una circunferencia: el unico foco es el centro")
            return [S._afin.centro_conica(c)]
        paso("Calculamos las tangentes desde I = ", self._i, ", y J = ", self._j)
        tangi = c.tangentes(self._i)
        _no_pasos()
        tangj = c.tangentes(self._j)
        _no_pasos(False)
        paso("Las de J deberian ser las conjugadas; desde I: ", tangi.ecuacion(factor = True), "; desde J: ", tangj.ecuacion(factor = True))
        (ri, si) = tangi.factorizacion()
        (rj, sj) = tangj.factorizacion()
        # Es parábola
        if ri == self.recta_infinito() or si == self.recta_infinito():
            r = ri if ri != self.recta_infinito() else si
            s = rj if rj != self.recta_infinito() else sj
            paso("Es una parabola: el foco sera la interseccion entre las rectas que no coinciden:")
            paso(r.implicitas()[0], ", ", s.implicitas()[0])
            _no_pasos()
            f = r.interseccion(s).punto().simplify_full()
            _no_pasos(False)
            return [f]
        else:
            _no_pasos()
            # Simplificamos por intentar quitar alguna i que sobre
            f1 = ri.interseccion(rj).punto().simplify_full()
            f2 = ri.interseccion(sj).punto().simplify_full()
            f3 = si.interseccion(rj).punto().simplify_full()
            f4 = si.interseccion(sj).punto().simplify_full()
            _no_pasos(False)
            paso("Intersecamos las cuatro rectas entre si:")
            paso(ri.implicitas()[0], " con: ", rj.implicitas()[0], ", ", sj.implicitas()[0], ": ", f1, ", ", f2)
            paso(si.implicitas()[0], ", con las mismas: ", f3, ", ", f4)
            # Devolvemos los focos reales primero
            return sorted([f1, f2, f3, f4], key = lambda f: not es_real(f))

    #\m
    # Devuelve una tupla con los ejes de la cónica dada. Para una cónica con cuatro focos, simplemente une los focos reales entre sí (y lo
    # mismo para los imaginarios). Para las cónicas con un solo foco, lo une con el "centro" infinito en el caso de la parábola y con el x
    # dado por parámetro en el caso de la circunferencia, pues el foco coincide con el centro (el otro eje es su perpendicular).
    #
    # IMPORTANTE. En el caso de la parábola el segundo eje que se devuelve no es un eje real, sino la recta del infinito. \\
    # IMPORTANTE 2. No es fiable 100% cuando depende de la factorización de cónicas degeneradas, pues no es un método fiable 100%.
    # Si se obtienen resultados muy raros o errores es recomendable ejecutar varias veces. Si eso no lo arregla, tocará calcular los focos a mano. \\
    # NOTA. No se muestra el procedimiento de calcular los focos.
    #
    # Parámetros \\
    # c: conica - cónica de la que obtener sus ejes \\
    # x: vector(3) - dirección (punto de la recta del infinito) desde la que calcular el prumer eje en el caso de la circunferencia
    # (por defecto ninguno, pues en general no es necesario)
    #
    def ejes(self, c, x = None):
        _no_pasos()
        focos = self.focos(c)
        _no_pasos(False)
        paso("Los focos son: ", focos)
        if len(focos) == 4:
            paso("Los unimos:")
            # Deberían estar ordenados, primero los reales
            e1 = subespacio(focos[0], focos[1])
            e2 = subespacio(focos[2], focos[3])
            return (e1, e2)
        elif not self.es_circunferencia(c):
            _no_pasos()
            centro = self._afin.puntos_infinitos_conica(c)[0]
            _no_pasos(False)
            paso("Lo unimos con el punto de tangencia en el infinito: ", centro)
            e1 = subespacio(focos[0], centro)
            paso("El otro sera la recta del infinito")
            return (e1, self.recta_infinito())
        else:
            assert x is not None and x in self.recta_infinito(), \
                "Para las circunferencias hay que dar un punto del infinito desde el que empezar"
            paso("Uno de los ejes es el que une: ", x, " con: ", focos[0], " (que es el centro)")
            e1 = subespacio(x, focos[0])
            paso("El otro eje es su perpendicular que pasa por: ", focos[0])
            e2 = self.perpendicular(self._afin.direccion(e1), focos[0])
            return (e1, e2)
        assert False, "Algo fue mal calculando los ejes =("

    #\m
    # Devuelve una lista con las directrices de la cónica indicada, esto es, las polares de los focos. Deberían ser rectas reales
    # o, en caso de ser imaginarias, simplificables. En el caso de que la cónica sea una circunferencia, su directriz será la propia
    # recta del infinito, luego no una recta del plano. Igualmente se devuelve el resultado.
    #
    # IMPORTANTE. No es fiable 100% cuando depende de la factorización de cónicas degeneradas, pues no es un método fiable 100%.
    # Si se obtienen resultados muy raros o errores es recomendable ejecutar varias veces. Si eso no lo arregla, tocará calcular los
    # focos a mano.
    #
    # Parámetros \\
    # c: conica - cónica de la que calcular sus directrices
    #
    def directrices(self, c):
        _no_pasos()
        if self.es_circunferencia(c):
            paso("Es una circunferencia, el resultado sera la recta del infinito, con lo que no es una directriz 'verdadera'")
        focos = self.focos(c)
        _no_pasos(False)
        paso("Los focos son: ", focos, "; calculamos sus polares")
        return map(lambda f: c.polar(f), focos)

    def __repr__(self):
        return "<Espacio euclideo de dimension 2 con recta del infinito de ecuacion " + str(self.recta_infinito().implicitas()[0]) + \
            " y puntos ciclicos conjugados del infinito I = " + str(self._i) + ", J = " + str(self._j) + ">"

#\f
# Determina si un punto es real.
#
# Implementación \\
# Comrueba que coincida con su conjugado.
#
# Parámetros \\
# p: vector(n) - punto a comprobar
#
def es_real(p):
    return bool(p == conjugate(p))
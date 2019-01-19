#\s
# Archivo con clases representantes de espacios tanto afines como euclídeos con diversas funciones relacionadas con
# aplicaciones, cónicas, etc.
#
# Autor: Pablo Sanz Sanz
#

# Clases

#\c
# Clase que representa un espacio afín de una dimensión arbitraria. Sirve para diversos cálculos relacionados con
# aplicaciones, cónicas, etc.
#
# NOTA. Se puede indicar que por defecto los resultados sean devueltos en la referencia canónica o no. En cualquier caso los
# puntos que se pasan como parámetro deben estar en la referencia de este espacio.
# Si, por ejemplo, queremos calcular el punto medio de los puntos (-1 : 0 : 1) y (1 : 0 : 1) (dados en la referencia canónica),
# pero la referencia de este espacio no es la canónica, hay que cambiar primero esos puntos a la referencia de este espacio a
# mano con la inversa de la matriz que devuelva el método matriz_cambio (M^-1 * X) y después llamar al método con estos
# dos nuevos puntos. Finalmente, en función del valor que se haya indicado a cambio_defecto, se cambiará la referencia del
# resultado no.
#
# Todo esto se aplica en general, salvo que se indique lo contrario.
#
class espacio_afin:

    #\i
    # Construye un espacio afín, dados el hiperplano del infinito y la referencia que se está utilizando.
    # Nótese que el hiperplano se supone que viene expresado en la referencia dada.
    #
    # Parámetros \\
    # hiperplano_infinito: subespacio(dim=n-1) - hiperplano que se está considerando como infinito (xn = 0 por defecto) \\
    # referencia: lista(vector(n+1))(long=n+2) - referencia en la que estarán expresados los puntos que se consideren en los
    # métodos de esta clase (canónica ([(1 : 0 : ... : 0), (0 : ... 0 : 1), (1 : ... : 1)]) por defecto) \\
    # cambio: booleano - indica si los valores devueltos deben ser previamente cambiados de referencia o no (True por defecto) \\
    # dimension: natural - dimensión de este espacio si no se han especificado los dos anteriores parámetros; por ejemplo,
    # espacio_afin(dimension = 2) crea un espacio afín cuya recta del infinto es z = 0 y la referencia es la canónica
    # (2 por defecto)
    #
    def __init__(self, hiperplano_infinito = None, referencia = None, cambio = True, dimension = 2):
        self._dim = hiperplano_infinito.dimension_ambiente() if hiperplano_infinito is not None else \
                    len(referencia[0]) - 1 if referencia is not None else \
                    dimension
        _no_pasos()
        # xn = 0 por defecto
        if hiperplano_infinito is None:
            hiperplano_infinito = subespacio(vector(([0] * (self._dim + 1)) + [1])).dual()
        _no_pasos(False)
        # I por defecto
        if referencia is None:
            referencia = referencia_canonica(self._dim).columns()
        assert hiperplano_infinito.dimension_ambiente() == self._dim, "La dimension ambiente del hiperplano debe coincidir con la de este espacio"
        assert hiperplano_infinito.dim() == self._dim - 1, "El subespacio debe ser un hiperplano"
        assert self.es_referencia_valida(referencia), "La referencia dada debe ser una referencia valida para este espacio afin"
        self._infinito = hiperplano_infinito
        self._referencia = referencia
        self._cambio = cambio
        m = matrix(referencia).T
        if m != 1:
            paso("Calculamos la matriz del cambio a la referencia canonica")
            self._matriz_cambio = matriz_asociada(m)
        else:
            _no_pasos()
            self._matriz_cambio = matriz_asociada(m)
            _no_pasos(False)

    #\m
    # Hace que los métodos cambien la referencia a la canónica por defecto. Si se usa s.cambio_defecto(True) (donde s es este espacio),
    # a partir de entonces las funciones devolverán los resultados cambiando primero la referencia.
    #
    # Parámetros \\
    # b: booleano - indica si se cambiará la referencia o no
    #
    def cambio_defecto(self, b):
        self._cambio = b

    # Métodos accedentes

    #\m
    # Devuelve la dimensión de este espacio afín.
    def dim(self):
        return self._dim

    #\m
    # Devuelve el hiperplano del infinito de este espacio afín en la referencia en que está expresado el espacio.
    def hiperplano_infinito(self):
        return self._infinito

    #\m
    # Devuelve la referencia en que está expresado este espacio.
    def referencia(self):
        return self._referencia

    #\m
    # Devuelve la matriz del cambio de referencia de la referencia en que se encuentra este espacio a la canónica.
    # Calcular la inversa para el cambio contrario.
    def matriz_cambio(self):
        return self._matriz_cambio

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
    # Determina si la referencia indicada (la propia de este espacio si no se indica ninguna) es una referencia afín, esto
    # es, cada punto de los n primeros pertenece al hiperplano del infinito.
    #
    # Parámetros \\
    # referencia: lista(vector(n+1))(long=n+2) - referencia a comprobar
    #
    def es_referencia_afin(self, referencia = None):
        if referencia is None:
            referencia = self._referencia
        else:
            assert self.es_referencia_valida(referencia), "La referencia dada debe ser valida"
        for p in referencia[0 : -2]:
            if p not in self._infinito:
                return False
        return True

    #\m
    # Devuelve un nuevo espacio afín (no modifica este) cuya referencia es la nueva indicada con el hiperplano del infinito
    # transformado correctamente. La variable que indica si se debe cambiar de referencia mantendrá su valor.
    #
    # NOTA. Se supone que los puntos dados están dados en coordenadas canónicas.
    #
    # Parámetros \\
    # referencia: lista(vector(n+1))(long=n+2) - nueva referencia del espacio
    #
    def cambiar_referencia(self, referencia):
        assert self.es_referencia_valida(referencia), "La referencia dada debe ser valida"
        ini = matrix(self._referencia).T
        fin = matrix(referencia).T
        paso("Obtenemos la matriz del cambio entre las referencias: ", self._referencia, ", ", referencia)
        m = cambiar_referencia(ini, fin)
        paso(m)
        paso("Cambiamos el hiperplano del infinito de referencia, cambiando de referencia el punto dual que representa con " + \
                "la matriz del cambio inversa y traspuesta M^-t = ", m.T^-1, " (resultado de Algebra)")
        _no_pasos()
        dual = self._infinito.dual().punto()
        coord = m.T^-1 * dual
        h = subespacio(coord).dual()
        _no_pasos(False)
        paso("Sus coordenadas: ", dual, " se convierten en: ", coord, ": ", h.implicitas()[0])
        return espacio_afin(h, referencia, self._cambio)

    #\m
    # Determina los puntos del infinito del subespacio dado, esto es la intersección de este con el hiperplano del infinito.
    # Devuelve un objeto del tipo subespacio. Si el subespacio dado es una recta se puede tomar el único punto de intersección
    # con el método punto() de la clase subespacio.
    #
    # NOTA. ¡No se cambia la referencia del subespacio devuelto! Habría que hacerlo a mano.
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
        return self.__p(r[-1])

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
        return self.__p(conjugado_armonico(p1, p2, pinf))

    # Métodos auxiliares

    # Para cambiar de referencia
    def __p(self, p):
        if self._cambio and self._matriz_cambio != 1:
            paso("Devolvemos el punto a la referencia canonica")
            return self._matriz_cambio * p
        return p

    def __repr__(self):
        return "<Espacio afin de dimension " + str(self._dim) + " con hiperplano del infinito de ecuacion " + \
                str(self._infinito.implicitas()[0]) + " y referencia " + str(self._referencia) + ">"
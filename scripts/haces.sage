#\s
# Archivo que contiene clases relativas a haces de hiperplanos y cónicas dados dos elementos del haz.
# Para usar haces de hiperplanos dada su base mejor usar subespacios duales.
#
# Autor: Pablo Sanz Sanz
#

# Clases


#\c
# Clase que representa un haz de hiperplanos (puntutales o duales).
#
class haz_hiperplanos:

    #\i
    # Crea un haz de hiperplanos dados dos hiperplanos representantes del haz. Será h1 + lambda * h2.
    #
    # Parámetros \\
    # h1: subespacio - primer hiperplano representante del haz \\
    # h2: subespacio - segundo hiperplano representante del haz
    def __init__(self, h1, h2):
        d = h1.dimension_ambiente()
        assert d == h2.dimension_ambiente(), "Los hiperplanos deben pertencer al mismo espacio"
        assert d - h1.dim() == 1 and d - h2.dim() == 1, "Los hiperplanos tienen codimension 1"
        self._param = var('lambda0', latex_name = '\\lambda')
        _no_pasos()
        self._generico = h1.dual().representantes()[0] + lambda0 * h2.dual().representantes()[0]
        _no_pasos(False)
        self._cero = h1
        self._infinito = h2

    # Métodos accedentes

    #\m
    # Devuelve la ecuación de este haz de hiperplanos. Nótese que puede quedar un denominador común, que se puede eliminar a mano si hiciera falta.
    def ecuacion(self):
        _no_pasos()
        res = subespacio(self._generico).dual().implicitas()[0]
        _no_pasos(False)
        return res

    #\m
    # Devuelve la dimensión del espacio en que se encuentra este haz.
    def dimension_ambiente(self):
        return self._cero.dimension_ambiente()

    #\m
    # Devuelve el hiperplano para la coordenada dada.
    #
    # Uso: H[lambda] (lambda es la coordenada del hiperplano y H un haz de hiperplanos).
    #
    # Si se quiere un hiperplano genérico del haz, sustituir por una variable (ej: r[x] devuelve un hiperplano dependiente de x).
    #
    # Parámetros \\
    # x: complejo/Infinity - coordenada del hiperplano que se quiere obtener
    #
    def __getitem__(self, x):
        if es_infinito(x):
            return self._infinito
        _no_pasos()
        res = subespacio(self._generico.substitute(self._param == x)).dual()
        _no_pasos(False)
        return res

    # Otros métodos

    #\m
    # Devuelve el subespacio base del haz, esto es, el subespacio intersección de los dos dados.
    #
    def base(self):
        paso("Intersecamos los dos hiperplanos dados al crear este haz")
        return self._cero.interseccion(self._infinito)

    #\m
    # Devuelve el hiperplano del haz que pasa por el punto dado.
    #
    # Implementación \\
    # Sustituye el punto en el subespacio genérico dependiente de lambda e iguala a 0.
    #
    # Parámetros \\
    # p: vector(n) - punto por el que pasa el hiperplano del haz
    #
    def forzar_punto(self, p):
        assert len(p) == self.dimension_ambiente() + 1, "El punto debe pertenecer al espacio del haz"
        _no_pasos()
        base = self.base()
        _no_pasos(False)
        assert p not in base, "No se puede forzar que el hiperplano pase por un punto base, pues todos pasan por ellos"
        ec = self._generico * p
        paso("Sustituimos el punto: ", p, " en: ", self.ecuacion())
        paso(ec, " = 0")
        # Resolvemos y nos quedamos con el resultado para sustituir
        sol = solve(ec == 0, self._param)
        if len(sol) == 0:
            return self[Infinity]
        return self[sol[0].rhs()]

    def __repr__(self):
        return "<Haz de hiperplanos de ecuacion " + str(self.ecuacion()) + ">"

#\c
# Clase que representa un haz de cónicas (puntutales o duales).
#
class haz_conicas:

    #\i
    # Crea un haz de cónicas dadas dos cónicas representantes del haz. Será c1 + lambda * c2.
    #
    # Parámetros \\
    # c1: cónica - primera cónica representante del haz \\
    # c2: cónica - segunda cónica representante del haz
    def __init__(self, c1, c2):
        self._param = var('lambda0', latex_name = '\\lambda')
        self._matriz = c1.matriz_asociada() + self._param * c2.matriz_asociada()
        self._vars = vector(var('x y z'))
        self._ecuacion = (self._vars * self._matriz * self._vars).factor()
        self._infinito = c2

    # Métodos accedentes

    #\m
    # Devuelve la matriz asociada a este haz de cónicas.
    def matriz_asociada(self):
        return self._matriz

    #\m
    # Devuelve la ecuación de este haz de cónicas.
    def ecuacion(self):
        return self._ecuacion == 0

    #\m
    # Devuelve la cónica para la coordenada dada.
    #
    # Uso: H[lambda] (lambda es la coordenada de la cónica y H un haz de cónicas).
    #
    # Si se quiere una cónica genérica del haz, sustituir por una variable (ej: r[x] devuelve una cónica dependiente de x).
    #
    # Parámetros \\
    # x: complejo/Infinity - coordenada de la cónica que se quiere obtener
    #
    def __getitem__(self, x):
        if es_infinito(x):
            return self._infinito
        return conica(self._matriz.substitute(self._param == x))

    # Otros métodos

    #\m
    # NO IMPLEMENTADO AÚN.
    #
    # Devuelve los puntos base del haz, esto es, los puntos por los que pasan todas las cónicas.
    #
    def base(self):
        return []

    #\m
    # Devuelve la cónica del haz que pasa por el punto dado.
    #
    # Implementación \\
    # Resuelve la ecuación p^t * H(lambda) * p = 0, donde H es la matriz representante del haz.
    #
    # Parámetros \\
    # p: vector(3) - punto por el que pasa la cónica del haz
    #
    def forzar_punto(self, p):
        _no_pasos()
        base = self.base()
        _no_pasos(False)
        assert p not in base, "No se puede forzar que la conica pase por un punto base, pues todas pasan por ellos"
        ec = p * self._matriz * p
        paso("Resolvemos:", matrix([p]), self._matriz, matrix([p]).T, " = ", ec, " = 0")
        # Resolvemos y nos quedamos con el resultado para sustituir
        sol = solve(ec == 0, self._param)
        if len(sol) == 0:
            return self[Infinity]
        return self[sol[0].rhs()]

    def __repr__(self):
        return "<Haz de conicas de ecuacion " + str(self.ecuacion()) + ">"

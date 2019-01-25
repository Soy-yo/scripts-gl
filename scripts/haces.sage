#\s
# Archivo que contiene clases relativas a haces de hiperplanos y cónicas dados dos elementos del haz.
# Para usar haces de hiperplanos dada su base mejor usar subespacios duales.
#
# Autor: Pablo Sanz Sanz
#

# Funciones globales

#\f
# Crea un haz genérico de cónicas, dados los cuatro puntos de contacto.
#
# Implementación \\
# Crea un haz con las cónicas degeneradas AB·CD y AC·BC.
#
# Parámetros \\
# a: vector(3) - primer punto de contacto \\
# b: vector(3) - segundo punto de contacto \\
# c: vector(3) - tercer punto de contacto \\
# d: vector(3) - cuarto punto de contacto
#
def haz_generico(a, b, c, d):
    assert len(a) == 3 and len(b) == 3 and len(c) == 3 and len(d) == 3, "Los puntos deben pertenecer al plano"
    _no_pasos()
    ab = subespacio(a, b)
    cd = subespacio(c, d)
    abcd = conica_degenerada(ab, cd)
    ac = subespacio(a, c)
    bd = subespacio(b, d)
    acbd = conica_degenerada(ac, bd)
    _no_pasos(False)
    paso("Creamos las conicas degeneradas AB*CD y AC*BD y generamos el haz:")
    paso("(", ab.implicitas()[0].lhs(), ")(", cd.implicitas()[0].lhs(), ") + lambda * (", ac.implicitas()[0].lhs(), ")(", bd.implicitas()[0].lhs(), ") = 0")
    return haz_conicas(abcd, acbd)

#\f
# Crea un haz con contacto simple, dados una cónica, dos puntos por donde pasa y el punto de contacto doble.
#
# Implementación \\
# Si A y B pertenecen a la cónica obtiene las rectas AT y BT, que formarán una cónica degenerada del haz. Con esta y la dada se construye el haz.
# SI no, se obtendrá también la tangente t a la cónica desde T de forma que t·AB es la otra cónica degenerada del haz.
#
# Parámetros \\
# conica: conica - cónica no degenerada que al menos contenga el punto T \\
# a: vector(3) - primer punto de contacto simple \\
# b: vector(3) - segundo punto de contacto simple \\
# t: vector(3) - punto de contacto doble
#
def haz_contacto_simple(conica, a, b, t):
    assert len(a) == 3 and len(b) == 3 and len(t) == 3, "Los puntos deben pertenecer al plano"
    assert t in conica, "El punto doble debe pertenecer a la conica"
    _no_pasos()
    at = subespacio(a, t)
    bt = subespacio(b, t)
    atbt = conica_degenerada(at, bt)
    _no_pasos(False)
    # Necesitamos también la tangente
    if a not in conica or b not in conica:
        paso("Calculamos la tangente desde T = ", t, " usando la polar")
        tang = conica.polar(t)
        _no_pasos()
        ab = subespacio(a, b)
        abt = conica_degenerada(tang, ab)
        _no_pasos(False)
        paso("Creamos las conicas degeneradas AT*BT y t*AB y generamos el haz:")
        paso("(", at.implicitas()[0].lhs(), ")(", bt.implicitas()[0].lhs(), ") + lambda * (", \
             tang.implicitas()[0].lhs(), ")(", ab.implicitas()[0].lhs(), ") = 0")
        return haz_conicas(atbt, abt)
    else:
        paso("Creamos la conica degenerada AT*BT y generamos el haz (xi es x, y, z respectivamente):")
        paso("(", conica.ecuacion().lhs(), ") + lambda * (", at.implicitas()[0].lhs(), ")(", bt.implicitas()[0].lhs(), ") = 0")
        return haz_conicas(conica, atbt)

#\f
# Crea un haz con contacto triple, dados una cónica, un punto de por el que pasa y el punto de contacto triple.
#
# NOTA. He visto que en algunos ejercicios se pide calcular el haz de contacto triple dado el punto A que no pertenece a la cónica.
# Como de momento no sé hacer eso, aquí no está permitido.
#
# Implementación \\
# Obtiene las rectas AT y t (tangente a la cónica desde T), que formarán una cónica degenerada del haz. Con esta y la dada se construye el haz.
#
# Parámetros \\
# conica: conica - una cónica del haz \\
# a: vector(3) - punto de contacto simple \\
# t: vector(3) - punto de contacto triple
#
def haz_contacto_triple(conica, a, t):
    assert len(a) == 3 and len(t) == 3, "Los puntos deben pertenecer al plano"
    assert a in conica and t in conica, "Los puntos deben pertenecer a la conica"
    paso("Calculamos la tangente desde T = ", t, " usando la polar")
    tang = conica.polar(t)
    _no_pasos()
    at = subespacio(a, t)
    att = conica_degenerada(at, tang)
    _no_pasos(False)
    paso("Creamos la conica degenerada AT*t y generamos el haz (xi es x, y, z respectivamente):")
    paso("(", conica.ecuacion().lhs(), ") + lambda * (", at.implicitas()[0].lhs(), ")(", tang.implicitas()[0].lhs(), ") = 0")
    return haz_conicas(conica, att)

#\f
# Crea un haz con contacto cuádruple, dados la cónica y el punto de contacto cuádruple.
#
# Implementación \\
# Obtiene la recta doble tangente a la cónica desde T, que formará una cónica degenerada del haz. Con esta y la dada se construye el haz.
#
# Parámetros \\
# conica: conica - una cónica no degenerada del haz \\
# t: vector(3) - punto de contacto cuádruple
#
def haz_contacto_cuadruple(conica, t):
    assert len(t) == 3, "El punto debe pertenecer al plano"
    assert t in conica, "El punto debe pertenecer a la conica"
    paso("Calculamos la tangente desde T = ", t, " usando la polar")
    tang = conica.polar(t)
    tt = conica_degenerada(tang)
    paso("Creamos la conica degenerada t^2 y generamos el haz (xi es x, y, z respectivamente):")
    paso("(", conica.ecuacion().lhs(), ") + lambda * (", tang.implicitas()[0].lhs(), ")^2 = 0")
    return haz_conicas(conica, tt)

#\f
# Crea un haz con dos contactos dobles, dados una o dos cónicas y los puntos de contacto doble.
#
# Implementación \\
# Si una cónica contiene los dos puntos S y T obtiene la recta doble ST, que formará una cónica degenerada del haz y
# se construye con esta y la dada. Si no se deben dar dos cónicas, una conteniendo S y la otra T, por lo que no
# pertenecerán al haz. Entonces se tomarán además las tangentes desde los puntos respectivos como segunda cónica
# degenerada del haz.
#
# Parámetros \\
# conica1: conica - una cónica del haz que contenga al menos a S o a T \\
# s: vector(3) - primer punto de contacto doble \\
# t: vector(3) - segundo punto de contacto doble \\
# conica2: conica - otra cónica del haz que contenga el punto de S y T que no contenga la primera, en caso de necesitarse
# (por defecto ninguna)
#
def haz_dos_contactos_dobles(conica1, s, t, conica2 = None):
    assert len(s) == 3 and len(t) == 3, "Los puntos deben pertenecer al plano"
    assert t in conica1 or s in conica1, "Alguno de los puntos debe pertenecer a la primera conica"
    assert t in conica1 and s in conica1 or conica2 is not None, "Los dos puntos deben estar en alguna conica del haz"
    _no_pasos()
    st = subespacio(s, t)
    sstt = conica_degenerada(st)
    _no_pasos(False)
    # Necesitamos también las tangentes
    if s not in conica1 or t not in conica1:
        assert s in conica2 or t in conica2, "Alguno de los puntos debe pertenecer a la segunda conica"
        # S será el punto de tangencia de la primera
        if s not in conica1:
            (s, t) = (t, s)
        paso("Calculamos la tangente desde S = ", s, " y T = ", t, " usando las polares")
        tang1 = conica1.polar(s)
        tang2 = conica2.polar(t)
        conica_tang = conica_degenerada(tang1, tang2)
        paso("Creamos las conicas degeneradas ST^2 y Tang(T)Tang(S) y generamos el haz:")
        paso("(", st.implicitas()[0].lhs(), ")^2 + lambda * (", tang1.implicitas()[0].lhs(), ")(", tang2.implicitas()[0].lhs(), ") = 0")
        return haz_conicas(sstt, conica_tang)
    else:
        paso("Creamos la conica degenerada ST^2 y generamos el haz (xi es x, y, z respectivamente):")
        paso("(", conica1.ecuacion().lhs(), ") + lambda * (", st.implicitas()[0].lhs(), ")^2 = 0")
        return haz_conicas(conica1, sstt)

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
    #
    def __init__(self, h1, h2):
        d = h1.dimension_ambiente()
        assert d == h2.dimension_ambiente(), "Los hiperplanos deben pertencer al mismo espacio"
        assert d - h1.dim() == 1 and d - h2.dim() == 1, "Los hiperplanos tienen codimension 1"
        self._param = var('lambda_var', latex_name = '\\lambda')
        _no_pasos()
        self._generico = h1.dual().punto() + self._param * h2.dual().punto()
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
    # Devuelve la variable simbólica que utilzia este haz como parámetro.
    def parametro(self):
        return self._param

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
# NOTA. No se implementa el método para obtener la base porque requeriría que una de las cónicas estuviera
# parametrizada. Si se quiere obtener la base deberá hacerse la intersección de dos cónicas representantes del haz
# a mano, eligiendo la parametrización de una de ellas.
#
class haz_conicas:

    #\i
    # Crea un haz de cónicas dadas dos cónicas representantes del haz. Será c1 + lambda * c2.
    #
    # Parámetros \\
    # c1: cónica - primera cónica representante del haz \\
    # c2: cónica - segunda cónica representante del haz
    def __init__(self, c1, c2):
        self._param = var('lambda_var', latex_name = '\\lambda')
        self._matriz = c1.matriz_asociada() + self._param * c2.matriz_asociada()
        variables = (var('x_var', latex_name = 'x'), var('y_var', latex_name = 'y'), var('z_var', latex_name = 'z'))
        self._vars = vector(variables)
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
    # Devuelve el vector de variables simbólicas que utiliza este haz para su representación.
    def variables(self):
        return self._vars

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

    #
    # Devuelve una tupla conteniendo los cuatro puntos base del haz, esto es, los puntos por los que pasan todas las cónicas.
    #
    # Implementación \\
    # Se obtiene una parametrización de una de las cónicas y se interseca con la otra directamente (ver conica.interseccion_conica
    # (cuadricas.sage)).
    #
    #def base(self):
    #    return []

    #\m
    # Devuelve la cónica del haz que pasa por el punto dado.
    #
    # NOTA. No se comprueba que el punto no esté en la base, pero se debe evitar para no obtener resultados inesperados.
    #
    # Implementación \\
    # Resuelve la ecuación p^t * H(lambda) * p = 0, donde H es la matriz representante del haz.
    #
    # Parámetros \\
    # p: vector(3) - punto por el que pasa la cónica del haz
    #
    def forzar_punto(self, p):
        #_no_pasos()
        #base = self.base()
        #_no_pasos(False)
        #assert all([matrix([p, x]).rank() == 2 for x in base]), "No se puede forzar que la conica pase por un punto base, pues todas pasan por ellos"
        ec = p * self._matriz * p
        paso("Resolvemos:", matrix([p]), self._matriz, matrix([p]).T, " = ", ec, " = 0")
        # Resolvemos y nos quedamos con el resultado para sustituir
        sol = solve(ec == 0, self._param)
        if len(sol) == 0:
            return self[Infinity]
        return self[sol[0].rhs()]

    #\m
    # Devuelve una lista conteniendo las hasta tres cónicas degeneradas de este haz.
    #
    # Implementación \\
    # Resuelve det(A0 + lambda*A1) = 0 para lambda, donde Ai son las matrices de las cónicas que generan el haz.
    # Después sustituye los resultados obtenidos para obtener las cónicas.
    #
    def degeneradas(self):
        ec = self._matriz.det().expand()
        paso("Resolvemos det", self._matriz, " = ", ec, " = 0")
        sol = [] if ec.degree(self._param) == 3 else [self._param == Infinity]
        sol = sol + solve(ec, self._param)
        paso(sol)
        paso("Obtenemos las conicas con esas coordenadas de: ", self._ecuacion)
        return map(lambda x: self[x.rhs()], sol)

    #\m
    # Devuelve la involución sobre una recta que genera este haz de cónicas, dada la recta.
    # No se comprueba, pero la recta no puede contener ningún punto base del haz.
    #
    # Implementación \\
    # Cada par de puntos en los que el haz corta la recta son los pares de la involución, con lo cual se interseca con las cónicas 0 e infinito y se
    # crea una homografía con esos dos pares, pues una involución queda determinada por dos pares.
    #
    # Parámetros \\
    # r: recta_proyectiva(dim=2) - recta sobre la que actuará la involución
    #
    def involucion_generada(self, r):
        assert r.dimension_ambiente() == 2, "La recta debe ser del plano"
        paso("Intersecamos la retca con las conicas de coordenada 0 e infinito")
        c1 = self[0]
        c2 = self[Infinity]
        paso("H[0] = ", c1.matriz_asociada(), "; H[", Infinity, "] = ", c2.matriz_asociada())
        (a, ap) = self[0].interseccion_recta(r)
        (b, bp) = self[Infinity].interseccion_recta(r)
        paso("Los pares son: ", (a, ap), " y ", (b, bp))
        return crear_homografia_recta(a, ap, b, bp, recta = r, involucion = True)

    #\m
    # Devuelve una tupla conteniendo las dos cónicas del haz tangentes a la recta dada.
    #
    # Implementación \\
    # Calcula los puntos fijos de la homografía que genera el haz con la recta dada y obtiene las cónicas que pasan por cada uno de ellos, pues cortan
    # la recta en puntos dobles.
    #
    # Parámetros \\
    # r: recta_proyectiva(dim=2) - recta de la que calcular las tangentes
    #
    def conicas_tangentes(self, r):
        assert r.dimension_ambiente() == 2, "La recta debe ser del plano"
        paso("Obtenemos la involucion que genera la recta:")
        h = self.involucion_generada(r)
        paso("Obtenemos los puntos fijos de la involucion")
        fijos = h.puntos_fijos()
        paso("Obtenemos las conicas que pasan por esos puntos fijos")
        return tuple(map(lambda x: self.forzar_punto(x), fijos))

    def __repr__(self):
        return "<Haz de conicas de ecuacion " + str(self.ecuacion()) + ">"

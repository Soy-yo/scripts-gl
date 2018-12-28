#\s
# Este archivo proporciona diferentes funciones y clases sobre razones dobles y homografías de la recta.
#
# Suponemos que los puntos expresados con coordenada inhomogénea theta se obtiene mediante x/y.
#
# Autor: Pablo Sanz Sanz
#

# Funciones globales

#\f
# Calcula la razón doble {A, B; C, D}. Es válido tanto para puntos como para parámetros inhomogéneos
# de una recta proyectiva (ver razon_doble_puntos y razon_doble_theta, respectivamente).
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto de la razón doble \\
# b: complejo/Infinity/vector(n>1) - segundo punto de la razón doble \\
# c: complejo/Infinity/vector(n>1) - tercer punto de la razón doble \\
# d: complejo/Infinity/vector(n>1) - cuarto punto de la razón doble
#
def razon_doble(a, b, c, d):
    if es_parametro(a):
        assert es_parametro(b) and es_parametro(c) and es_parametro(d), \
                "Todos los puntos deben ser del mismo tipo"
        return razon_doble_theta(a, b, c, d)
    return razon_doble_puntos(a, b, c, d)

#\f
# Calcula la razón doble ({P0, P1, P2, P3}) de los puntos alineados dados.
#
# Implementación \\
# Utiliza que la razón doble es invariante por proyecciones y la fórmula de la
# razón doble de los determinantes. Se elige una proyección que siga siendo una recta.
#
# Parámetros \\
# p0: vector(n>1) - primer punto de la razón doble \\
# p1: vector(n>1) - segundo punto de la razón doble \\
# p2: vector(n>1) - tercer punto de la razón doble \\
# p3: vector(n>1) - cuarto punto de la razón doble
#
def razon_doble_puntos(p0, p1, p2, p3):
    assert len(p0) > 1, "Para coordenadas de la recta usar 'razon_doble'"
    _no_pasos()
    assert subespacio(p0, p1, p2, p3).dim() == 1, "Los puntos deben estar alineados"
    _no_pasos(False)
    fin = len(p0) == 2
    i = 0
    j = 1
    q0 = p0
    q1 = p1
    q2 = p2
    q3 = p3
    while not fin:
        # Proyectamos en las coordenadas i y j
        q0 = vector([p0[i], p0[j]])
        q1 = vector([p1[i], p1[j]])
        q2 = vector([p2[i], p2[j]])
        q3 = vector([p3[i], p3[j]])
        # Si nos sigue quedando una recta (y no un punto) ya podemos calcular la razón doble
        _no_pasos()
        fin = q0 != 0 and q1 != 0 and q2 != 0 and q3 != 0 and subespacio(q0, q1, q2, q3).dim() == 1
        _no_pasos(False)
        j = j + 1
        # Creo que ni debería cumplirse esta condición
        if j == len(p0):
            i = i + 1
            j = i + 1
    if len(p0) > 2:
        paso("Como la razon doble es invariante por proyecciones podemos usar los siguientes vectores:")
        paso(q0, ", ", q1, ", ", q2, ", ", q3)
    m0 = matrix([[q0[0], q2[0]], [q0[1], q2[1]]])
    m1 = matrix([[q1[0], q2[0]], [q1[1], q2[1]]])
    m2 = matrix([[q0[0], q3[0]], [q0[1], q3[1]]])
    m3 = matrix([[q1[0], q3[0]], [q1[1], q3[1]]])
    det0 = m0.det()
    det1 = m1.det()
    det2 = m2.det()
    det3 = m3.det()
    paso("Calculamos los determinantes:")
    paso("det", m0, " = ", det0)
    paso("det", m1, " = ", det1)
    paso("det", m2, " = ", det2)
    paso("det", m3, " = ", det3)
    paso("Calculamos la razon doble: {P0, P1; P2, P3} = (", det0, "/", det1, ") : (", det2, "/", det3, ")")
    if det1 == 0 or det2 == 0:
        return Infinity
    return (det0 / det1) * (det3 / det2)

#\f
# Calcula la razón doble ({theta0, theta1, theta2, theta3}) de los puntos de la recta dados
# (con parámetro no homogéneo).
#
# Implementación \\
# Utiliza la razón doble definida para puntos tras homogeneizar las coordenadas.
#
# Parámetros \\
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble \\
# theta3: complejo/Infinity - cuarto punto de la razón doble
#
def razon_doble_theta(theta0, theta1, theta2, theta3):
    p0 = vector([1, 0]) if es_infinito(theta0) else vector([theta0, 1])
    p1 = vector([1, 0]) if es_infinito(theta1) else vector([theta1, 1])
    p2 = vector([1, 0]) if es_infinito(theta2) else vector([theta2, 1])
    p3 = vector([1, 0]) if es_infinito(theta3) else vector([theta3, 1])
    paso("Homogeneizando las coordenadas quedan los vectores:")
    paso(p0, ", ", p1, ", ", p2, ", ", p3)
    return razon_doble_puntos(p0, p1, p2, p3)

#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es, {A, B; C, D} = -1. Es válido tanto
# para puntos como para parámetros inhomogéneos de una recta proyectiva (ver es_cuaterna_armonica_puntos y
# es_cuaterna_armonica_theta, respectivamente).
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto de la razón doble \\
# b: complejo/Infinity/vector(n>1) - segundo punto de la razón doble \\
# c: complejo/Infinity/vector(n>1) - tercer punto de la razón doble \\
# d: complejo/Infinity/vector(n>1) - cuarto punto de la razón doble
#
def es_cuaterna_armonica(a, b, c, d):
    if es_parametro(a):
        assert es_parametro(b) and es_parametro(c) and es_parametro(d), \
                "Todos los puntos deben ser del mismo tipo"
        return es_cuaterna_armonica_theta(a, b, c, d)
    return es_cuaterna_armonica_puntos(a, b, c, d)


#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es, {P0, P1, P2, P3} = -1.
#
# Implementación \\
# Devuelve {P0, P1; P2, P3} == -1 utilizando las funciones anteriores.
#
# Parámetros \\
# p0: vector(n>1) - primer punto de la razón doble \\
# p1: vector(n>1) - segundo punto de la razón doble \\
# p2: vector(n>1) - tercer punto de la razón doble \\
# p3: vector(n>1) - cuarto punto de la razón doble
#
def es_cuaterna_armonica_puntos(p0, p1, p2, p3):
    return razon_doble_puntos(p0, p1, p2, p3) == -1

#\f
# Determina si los puntos dados forman una cuaterna armónica, esto es,
# {theta0, theta1, theta2, theta3} = -1.
#
# Implementación \\
# Devuelve {theta0, theta1; theta2, theta3} == -1 utilizando las funciones anteriores.
#
# Parámetros \\
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble \\
# theta3: complejo/Infinity - cuarto punto de la razón doble
#
def es_cuaterna_armonica_theta(theta0, theta1, theta2, theta3):
    return razon_doble(theta0, theta1, theta2, theta3) == -1

#\f
# Devuelve el conjugado armónico de c respecto de a y b. Es válido tanto para puntos como para parámetros inhomogéneos
# de una recta proyectiva (ver conjugado_armonico_puntos y conjugado_armonico_theta, respectivamente).
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto de la razón doble \\
# b: complejo/Infinity/vector(n>1) - segundo punto de la razón doble \\
# c: complejo/Infinity/vector(n>1) - tercer punto de la razón doble \\
#
def conjugado_armonico(a, b, c):
    if es_parametro(a):
        assert es_parametro(b) and es_parametro(c), "Todos los puntos deben ser del mismo tipo"
        return conjugado_armonico_theta(a, b, c)
    return conjugado_armonico_puntos(a, b, c)

#\f
# Devuelve el conjugado armónico de p2 respecto de p0 y p1.
#
# Implementación \\
# Como la razón doble es la coordenada del cuarto punto respecto a la referencia {p0, p1; p2}, obtiene la recta que
# forman p0, p1 y p2 y devuelve el punto correspondiente a la coordenada theta == -1.
#
# Parámetros \\
# p0: vector(n>1) - primer punto de la razón doble \\
# p1: vector(n>1) - segundo punto de la razón doble \\
# p2: vector(n>1) - tercer punto de la razón doble
#
def conjugado_armonico_puntos(p0, p1, p2):
    paso("Montamos la recta con la referencia {P0=", p0, ", P1=", p1, "; P2=", p2, "}")
    paso("Asi, P0: theta = Infinity, P1: theta = 0, P2: theta = 1 y sustituimos en -1")
    return recta_proyectiva(p0, p1, p2)[-1]

#\f
# Devuelve el conjugado armónico de theta2 respecto de theta0 y theta1.
#
# Implementación \\
# Resuelve la ecuación {theta0, theta1, theta2, theta} == -1 en theta utilizando las funciones anteriores.
#
# Parámetros \\
# theta0: complejo/Infinity - primer punto de la razón doble \\
# theta1: complejo/Infinity - segundo punto de la razón doble \\
# theta2: complejo/Infinity - tercer punto de la razón doble
#
def conjugado_armonico_theta(theta0, theta1, theta2):
    theta = var('theta', latex_name = r'theta')
    paso("Resolvemos {", theta0, ", ", theta1, ", ", theta2, ", ", theta, "} = -1")
    res = solve([es_cuaterna_armonica_theta(theta0, theta1, theta2, theta)], theta)
    # No ha habido solución
    if len(res) == 0:
        return Infinity
    return res[0].rhs()

#\f
# Crea un objeto del tipo homografia_recta dados la recta y tres puntos y sus imágenes. Es válido tanto para puntos
# como parámetros inhomogéneos. \\
# Si se quiere crear una involución dados dos puntos y sus imágenes se pueden omitir los terceros puntos. La función se podría
# llamar de la siguiente forma: h = crear_homografia_recta(a, ap, b, bp, recta = r, involucion = True), si se especifica una recta r.
#
# Diferenciar de crear_aplicacion_proyectiva (aplicaciones.sage) en el orden de los argumentos y en que aquí no se piden listas.
#
# Implementación \\
# Si se dan puntos se obtiene primero su coordenada de la recta. Una vez obtenidas las coordenadas se usa
# crear_homografia_recta_theta.
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto a transformar \\
# ap: complejo/Infinity/vector(n>1) - imagen del primer punto \\
# b: complejo/Infinity/vector(n>1) - segundo punto a transformar \\
# bp: complejo/Infinity/vector(n>1) - imagen del segundo punto \\
# c: complejo/Infinity/vector(n>1) - tercer punto a transformar (no necesario si se va a crear una involución) \\
# cp: complejo/Infinity/vector(n>1) - imagen del tercer punto (no necesario si se va a crear una involución) \\
# recta: recta_proyectiva - recta sobre la que actúa la homografía (por defecto la recta que ponga homografia_recta
# por defecto (sólo válido para parámetros inhomogéneos)) \\
# involucion: booleano - True si se quiere crear una involución, False si no necesariamente (False por defecto)
#
def crear_homografia_recta(a, ap, b, bp, c = None, cp = None, recta = None, involucion = False):
    assert c is not None and cp is not None or involucion, "Si no se crea una involucion se deben indicar los terceros puntos"
    if es_parametro(a):
        assert es_parametro(ap) and es_parametro(b) and es_parametro(bp), "Todos los puntos deben ser del mismo tipo"
        # Si no es una involución también hay que comprobar el tercer punto
        if not involucion:
            assert es_parametro(c) and es_parametro(cp), "Todos los puntos deben ser del mismo tipo"
        return crear_homografia_recta_theta(a, ap, b, bp, c, cp, recta, involucion)
    assert recta != None, "Si no se especifica la recta se deben usar parametros inhomogeneos"
    theta0 = recta.coordenada(a)
    theta0p = recta.coordenada(ap)
    theta1 = recta.coordenada(b)
    theta1p = recta.coordenada(bp)
    theta2 = recta.coordenada(c) if not involucion else None
    theta2p = recta.coordenada(cp) if not involucion else None
    paso("Las coordenadas en la recta de los puntos que se van a transformar son:")
    paso(theta0, " -> ", theta0p)
    paso(theta1, " -> ", theta1p)
    if not involucion:
        paso(theta2, " -> ", theta2p)
    return crear_homografia_recta_theta(theta0, theta0p, theta1, theta1p, theta2, theta2p, recta, involucion)

#\f
# Crea un objeto del tipo homografia_recta dados la recta y tres puntos (con parámetros inhomogéneos) y sus imágenes. \\
# Si se quiere crear una involución dados dos puntos y sus imágenes se pueden omitir los terceros puntos. La función se podría
# llamar de la siguiente forma: h = crear_homografia_recta_theta(theta0, theta0p, theta1, theta1p, recta = r, involucion = True),
# si se especifica una recta r.
#
# Diferenciar de crear_aplicacion_proyectiva (aplicaciones.sage) en el orden de los argumentos y en que aquí no se piden listas.
#
# Implementación \\
# Se crea una primera homografía con ecuación theta' = (a*theta + b) / (c*theta + d), se sustituyen cada uno de los
# puntos dados y se resuelve el sistema.
#
# Parámetros \\
# theta0: complejo/Infinity - primer punto a transformar \\
# theta0p: complejo/Infinity - imagen del primer punto \\
# theta1: complejo/Infinity - segundo punto a transformar \\
# theta1p: complejo/Infinity - imagen del segundo punto \\
# theta2: complejo/Infinity - tercer punto a transformar (no necesario si se va a crear una involución) \\
# theta2p: complejo/Infinity - imagen del tercer punto (no necesario si se va a crear una involución) \\
# recta: recta_proyectiva - recta sobre la que actúa la homografía (por defecto la recta que ponga homografia_recta
# por defecto (sólo válido para parámetros inhomogéneos)) \\
# involucion: booleano - True si se quiere crear una involución, False si no necesariamente (False por defecto)
#
def crear_homografia_recta_theta(theta0, theta0p, theta1, theta1p, theta2 = None, theta2p = None, recta = None, involucion = False):
    assert theta2 is not None and theta2p is not None or involucion, "Si no se crea una involucion se deben indicar los terceros puntos"
    var('a b c d')
    rhs0 = a / c if es_infinito(theta0) else (a * theta0 + b) / (c * theta0 + d)
    rhs1 = a / c if es_infinito(theta1) else (a * theta1 + b) / (c * theta1 + d)
    if not involucion:
        rhs2 = a / c if es_infinito(theta2) else (a * theta2 + b) / (c * theta2 + d)
    ec0 = rhs0.denominator() == 0 if es_infinito(theta0p) else theta0p == rhs0
    ec1 = rhs1.denominator() == 0 if es_infinito(theta1p) else theta1p == rhs1
    if not involucion:
        ec2 = rhs2.denominator() == 0 if es_infinito(theta2p) else theta2p == rhs2
    # Una involución tiene traza nula
    sistema = [ec0, ec1, ec2] if not involucion else [ec0, ec1, a == -d]
    sol = solve(sistema, a, b, c, d)[0]
    m11 = sol[0].rhs()
    m12 = sol[1].rhs()
    m21 = sol[2].rhs()
    m22 = sol[3].rhs()
    div = gcd([m11, m12, m21, m22])
    paso("La homografia tendra la forma: theta' = (a theta + b) / (c theta + d)")
    if involucion:
        paso("Para que la homografia sea una involucion debemos asegurar traza = 0, luego a = -d")
    paso("Resolvemos el sistema:", sistema, "para a, b, c, d:")
    paso([a == m11 / div, b == m12 / div, c == m21 / div, d == m22 / div])
    m = matrix([[m11 / div, m12 / div], [m21 / div, m22 / div]])
    _no_pasos()
    res = homografia_recta(m, recta)
    _no_pasos(False)
    return res

#\f
# Crea una homografía entre dos rectas que se cortan dados tres puntos de la recta origen y sus imágenes en la recta destino.
# Es válido tanto para puntos como parámetros inhomogéneos.
#
# Implementación \\
# Con las coordenadas inhomogéneas de cada punto se halla la homografía sobre una falsa recta {Infinity, 0; 1} que los
# trasnforma de la forma indicada. Simplemente se transforma la coordenada de una recta a la coordenada de la otra como si fuera
# una sola recta.
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - primer punto a transformar \\
# ap: complejo/Infinity/vector(n>1) - imagen del primer punto \\
# b: complejo/Infinity/vector(n>1) - segundo punto a transformar \\
# bp: complejo/Infinity/vector(n>1) - imagen del segundo punto \\
# c: complejo/Infinity/vector(n>1) - tercer punto a transformar \\
# cp: complejo/Infinity/vector(n>1) - imagen del tercer punto \\
# r_origen: recta_proyectiva - recta desde donde parte la homografía \\
# r_destino: recta_proyectiva - recta de llegada de la homografía
#
def crear_homografia_dos_rectas(a, ap, b, bp, c, cp, r_origen, r_destino):
    if es_parametro(a):
        assert es_parametro(ap) and es_parametro(b) and es_parametro(bp) and es_parametro(c) and es_parametro(cp), \
                "Todos los puntos deben ser del mismo tipo"
        paso("Creamos una homografia sobre una falsa recta que convierta: ", a, " -> ", ap, ",", b, " -> ", bp, ",", c, " -> ", cp)
        h = crear_homografia_recta_theta(a, ap, b, bp, c, cp)
        return homografia_dos_rectas(h, r_origen, r_destino)
    # No es un parámetro, será un punto
    assert a in r_origen and b in r_origen and c in r_origen, "Los puntos a, b y c deben pertenecer a la recta origen"
    assert ap in r_destino and bp in r_destino and cp in r_destino, "Los puntos a', b' y c' deben pertenecer a la recta destino"
    theta0 = r_origen.coordenada(a)
    theta0p = r_destino.coordenada(ap)
    theta1 = r_origen.coordenada(b)
    theta1p = r_destino.coordenada(bp)
    theta2 = r_origen.coordenada(c)
    theta2p = r_destino.coordenada(cp)
    paso("Creamos una homografia sobre una falsa recta que convierta: ", \
         theta0, " -> ", theta0p, ",", theta1, " -> ", theta1p, ",", theta2, " -> ", theta2p)
    h = crear_homografia_recta_theta(theta0, theta0p, theta1, theta1p, theta2, theta2p)
    return homografia_dos_rectas(h, r_origen, r_destino)

#\f
# Crea una homografía entre dos rectas que se cortan dados un punto y su imagen y el eje de la homografía (como subespacio).
# Es válido tanto para puntos como parámetros inhomogéneos.
#
# Implementación \\
# Obtiene el punto de corte del eje con cada una de las rectas. En caso de ser dos puntos distintos X e Y se tendrá que O' = Y
# y X' = O. Si son el mismo punto se utilizará un punto arbitrario extra B, que se puede especificar si se quiere, y se calculará
# su imagen: se traza una recta que pase por B y A' y cortará el eje en P; uniendo A con P se obtiene B' en la recta destino.
# El par de puntos será ese corte doble del eje con las rectas. Una vez obtenidos los puntos se usa crear_homografia_dos_rectas().
#
# NOTA. Por simplicidad no se hacen comprobaciones sobre puntos coincidentes, pero se debe asegurar que ninguno de los puntos
# que se da como parámetros (tanto a, ap, como extra) coincida con la intersección (o entre ellos). En caso contrario es probable
# encontrarse con resultados extraños o incluso errores.
#
# Parámetros \\
# a: complejo/Infinity/vector(n>1) - punto a transformar \\
# ap: complejo/Infinity/vector(n>1) - imagen del punto \\
# r_origen: recta_proyectiva - recta desde donde parte la homografía \\
# r_destino: recta_proyectiva - recta de llegada de la homografía \\
# eje: subespacio - eje de la homografía, que debe cortar a las anteriores rectas (nótese que se pide un subespacio) \\
# extra: complejo/Infinity/vector(n>1) - punto extra de r_origen que se utilizará en caso de que el corte del eje con las rectas
# sea un único punto (por defecto la suma del corte y el punto a dado, como puntos)
#
def crear_homografia_dos_rectas_eje(a, ap, r_origen, r_destino, eje, extra = None):
    if es_parametro(a):
        assert es_parametro(ap), "Todos los puntos deben ser del mismo tipo"
        a = r_origen[a]
        ap = r_destino[ap]
    _no_pasos()
    x = eje.interseccion(r_origen.subespacio())
    y = eje.interseccion(r_destino.subespacio())
    o = r_origen.subespacio().interseccion(r_destino.subespacio())
    _no_pasos(False)
    assert o.es_punto(), "La interseccion de las rectas debe ser un unico punto"
    # Distinto punto de corte
    if x != y:
        b = x.representantes()[0]
        bp = o.representantes()[0]
        c = bp
        cp = y.representantes()[0]
        paso("Los puntos de corte del eje con las rectas son, respectivamente, ", b, ", ", cp)
        return crear_homografia_dos_rectas(a, ap, b, bp, c, cp, r_origen, r_destino)
    # Mismo punto de corte
    else:
        b = x.representantes()[0]
        c = a + b if extra is None else extra
        # Comprobemos que el punto extra no coincide con ninguno de los otros
        paso("El punto de corte del eje con las rectas ha sido unico:", b)
        paso("Hallamos la imagen del punto extra escogido:", c, ", uniendolo con A':", ap, " y cortando esta recta con el eje")
        _no_pasos()
        r = subespacio(c, ap)
        z = r.interseccion(eje).representantes()[0]
        s = subespacio(a, z)
        cp = s.interseccion(r_destino.subespacio()).representantes()[0]
        _no_pasos(False)
        paso("Resumiendo, la recta corta en el eje en: ", z, ", y uniendolo con A corta en la recta destino en: ", cp)
        return crear_homografia_dos_rectas(a, ap, b, b, c, cp, r_origen, r_destino)

# Clases

#\c
# Clase que representa una recta proyectiva para cualquier dimension ambiente. Sirve para cálculos de razones dobles
# o para homografías entre rectas.
#
class recta_proyectiva:

    #\i
    # Inicializa la recta dada la referencia {p0, p1; p2} que se quiere utilizar.
    #
    # Implementación \\
    # Obtiene la matriz asociada a la referencia dada (esto es, multiplicando el vector asociado a p0 y p1 por un coeficiente
    # de forma que a*p0 + b*p1 == p2) y construye el punto b*p1 + theta*a*p0.
    #
    # Parámetros \\
    # p0: vector(n) - punto cuya coordenada será theta == Infinity \\
    # p1: vector(n) - punto cuya coordenada será theta == 0 \\
    # p2: vector(n) - punto cuya coordenada será theta == 1
    #
    def __init__(self, p0, p1, p2):
        _no_pasos()
        self._subespacio = subespacio(p0, p1, p2)
        _no_pasos(False)
        assert self._subespacio.dim() == 1, "Los puntos deben estar alineados"
        # Forzamos p2: theta == 1
        var('a b')
        paso("Forzamos P2: theta = 1:")
        paso(a, "*", p0, " + ", b, "*", p1, " = ", a * p0 + b * p1, " = ", p2)
        coef = solve((a * p0 + b * p1 - p2).list(), a, b)
        paso(coef[0])
        self._theta = var('theta', latex_name = r'theta')
        self._infinito = p0
        self._valores = (b * p1 + self._theta * a * p0).substitute(coef[0])
        paso("Sustituimos y combinamos con theta: ", self._valores)

    # Métodos accedentes

    #\m
    # Devuelve el punto con coordenada theta.
    #
    # Uso: r[theta_0] (theta_0 es la coordenada del punto y r es una recta(_proyectiva)).
    #
    # Parámetros \\
    # theta: complejo/Infinity - coordenada del punto de la recta que se quiere obtener
    #
    def __getitem__(self, theta):
        if es_infinito(theta):
            return self._infinito
        return self._valores.substitute(self._theta == theta)

    #\m
    # Devuelve la referencia de esta recta proyectiva en forma de lista [Infinity, 0, 1].
    def referencia(self):
        return [self[Infinity], self[0], self[1]]

    #\m
    # Devuelve esta recta proyectiva como objeto del tipo subespacio.
    def subespacio(self):
        return self._subespacio

    #\m
    # Devuelve la dimensión del espacio en que se encuentra esta recta.
    def dimension_ambiente(self):
        return self._subespacio.dimension_ambiente()

    # Otros métodos

    #\m
    # Devuelve la coordenada no homogénea asociada al punto p dado.
    #
    # Implementación \\
    # Resuelve la ecuación recta(theta) == alpha * p para theta.
    #
    # Parámetros \\
    # p: vector(n) - punto de la recta del que se quiere conocer su coordenada
    #
    def coordenada(self, p):
        assert p in self, "El punto debe pertenecer a la recta"
        var('alpha')
        inf = len(solve([self._infinito[i] == alpha * p[i] for i in range(len(p))], alpha)) > 0
        if inf:
            return Infinity
        return solve([self._valores[i] == alpha * p[i] for i in range(len(p))], alpha, self._theta)[0][1].rhs()

    # \m
    # Operador ==. Determina si dos rectas proyectivas son iguales.
    #
    # Uso: r == s (r y s son rectas proyectivas).
    #
    # Implementación \\
    # Comrpueba que el subespacio sea el mismo.
    #
    # Parámetros \\
    # otra: recta_proyectiva - recta a comprobar la igualdad
    #
    def __eq__(self, otra):
        return self._subespacio == otra._subespacio

    # \m
    # Operador in. Determina si un punto está contenido en este subespacio o no.
    #
    # Uso: P in r (P es un punto y r una recta proyectiva).
    #
    # Implementación \\
    # Comprueba que el punto pertenezca al subespacio asociado.
    #
    # Parámetros \\
    # punto: vector(n) - punto que comprobar si pertenece a la recta
    #
    def __contains__(self, punto):
        return punto in self._subespacio

    def __repr__(self):
        return "<Recta proyectiva " + str(self._valores) + ">"

#\c
# Clase que representa una homografía de una recta en sí misma. En general, se crearán usando funciones auxiliares
# y no con su constructor.
#
class homografia_recta:

    #\i
    # Construye una homografía sobre la recta dada cuya matriz sea la dada. En general, este constructor no se
    # usará directamente.
    #
    # Parámetros \\
    # matriz: matriz(2, 2) - matriz asociada a la homografía \\
    # recta: recta_proyectiva - recta sobre la que actúa la homografía (por defecto una recta con referencia {Infinity, 0; 1})
    #
    def __init__(self, matriz, recta = None):
        assert matriz.nrows() == 2 and matriz.ncols() == 2, "La matriz debe ser 2x2"
        assert matriz.det() != 0, "Para que sea una homografia el rango de su matriz asociada debe ser maximo"
        if recta is None:
            _no_pasos()
            recta = recta_proyectiva(vector([1, 0]), vector([0, 1]), vector([1, 1]))
            _no_pasos(False)
        self._matriz = matriz
        self._recta = recta
        self._theta1 = var('theta', latex_name = r'theta')
        self._theta2 = var('theta2', latex_name = r'theta_2')
        self._expresion_mobius = theta2 == (matriz[0][0] * theta + matriz[0][1]) / (matriz[1][0] * theta + matriz[1][1])
        # Multiplicamos por el denominador y restamos el numerador
        self._expresion = (theta2 * self._expresion_mobius.rhs().denominator() - self._expresion_mobius.rhs().numerator()).expand() == 0

    # Métodos accedentes

    #\m
    # Devuelve esta homografía como una aplicacion_proyectiva.
    def aplicacion(self):
        return aplicacion_proyectiva(self._matriz)

    #\m
    # Devuelve la matriz asociada a esta homografía.
    def matriz_asociada(self):
        return self._matriz

    #\m
    # Devuelve la recta sobre la que se aplica esta homografía.
    def recta(self):
        return self._recta

    #\m
    # Devuelve la referencia de la recta en la que está expresada esta homografía.
    def referencia(self):
        return self._recta.referencia()

    #\m
    # Devuelve la homografía expresada como una transformación de Möbius.
    def expresion_mobius(self):
        return self._expresion_mobius

    #\m
    # Devuelve la ecuación en theta y theta' de esta homografía.
    def ecuacion(self):
        return self._expresion

    # Otros métodos

    #\m
    # Calcula la imagen mediante esta homografía del punto dado.
    #
    # Uso: h(x) ó h(p) (donde r es una homografia_recta, x un complejo/Infinity y p un punto).
    #
    # Implementación \\
    # Dado un parámetro no homogéneo sustituye en la expresión de Möbius y devuelve el resultado. \\
    # Dado un punto del espacio, primero calcula su coordenada en la recta y luego vuelve a llamar a este método.
    #
    # Parámetros \\
    # x: complejo/Infinity/vector(n) - punto del que se quiere calcular su imagen
    #
    def __call__(self, x):
        # Estamos recibiendo un parámetro no homogéneo, simplemente sustituimos
        if es_parametro(x):
            return self.__sustituir(x)
        # Asumimos que es un vector, si es del mismo espacio calculamos su coordenada en la recta primero
        assert x in self._recta.subespacio(), "El punto debe pertenecer a la recta"
        coord = self._recta.coordenada(x)
        paso("La coordenada de ", x, " en la recta es: ", coord, ". Sustituimos en la expresion y recuperamos el punto")
        return self._recta[self(coord)]

    #\m
    # Devuelve los autovalores de esta homografía.
    #
    # Implementación \\
    # Ver autovalores() de aplicacion_proyectiva (aplicacion_proyectiva.sage).
    #
    def autovalores(self):
        return self.aplicacion().autovalores()

    #\m
    # Devuelve los parámetros de los puntos fijos de esta homografía. Devuelve una lista vacía en caso de que sea la identidad.
    #
    # Implementación \\
    # Calcula los autovectores y deshomogeniza.
    #
    def puntos_fijos_theta(self):
        paso("Podemos obtener los puntos fijos directamente de:", self._expresion.substitute(self._theta1 == self._theta2), \
                ", pero vamos a usar los autovectores porque es más sencillo para el algoritmo")
        if self.es_identidad():
            return []
        return map(lambda x: Infinity if x[1] == 0 else x[0] / x[1], self.aplicacion().puntos_fijos())

    #\m
    # Devuelve los puntos fijos de esta homografía. Devuelve una lista vacía en caso de que sea la identidad.
    #
    # Implementación \\
    # Calcula los autovectores de la matriz asociada.
    #
    def puntos_fijos(self):
        fijos = self.puntos_fijos_theta()
        paso("Sustituimos los puntos fijos ", fijos, "en la recta")
        return map(lambda theta: self._recta[theta], fijos)

    #\m
    # Determina si esta homografía es la identidad.
    #
    # Implementación \\
    # Sustituye theta' por theta y resta un lado de la expresión de Möbius al otro para comparar con 0.
    #
    def es_identidad(self):
        return bool(self._expresion_mobius.lhs().substitute(self._theta2 == self._theta1) - self._expresion_mobius.rhs() == 0)

    #\m
    # Determina si esta homografía es una elación, esto es si tiene un punto fijo doble.
    #
    # Implementación \\
    # Calcula los puntos fijos y compara.
    #
    def es_elacion(self):
        return len(self.puntos_fijos()) == 1

    #\m
    # Determina si esta homografía es una involución, esto es si self^2 == self.
    #
    # Implementación \\
    # Comprueba que la traza sea nula.
    #
    def es_involucion(self):
        return self._matriz.trace() == 0

    #\m
    # Calcula el módulo de esta homografía, esto es el cociente de sus autovalores (en cualquier orden).
    #
    def modulo(self):
        paso("El módulo se calcula mediante el cociente de los autovalores:")
        autovalores = self.autovalores()
        paso("Autovalores:", autovalores)
        # Autovalor doble: la división siempre da 1
        if len(autovalores) == 1:
            return 1
        return autovalores[0] / autovalores[1]

    #\m
    # Devuelve esta misma homografía en una referencia adecuada de forma que su ecuación quede de la forma más
    # simple posible. Para las involuciones no devuelve la forma general theta*theta' = 1. Para ello usar simplificar_involucion().
    #
    # Implementación \\
    # En función de los puntos fijos devuelve: \\
    # infinitos -> ella misma (identidad) \\
    # 1         -> matriz (1, 1; 0, 1) con referencia el punto fijo y dos arbitrarios, multiplicada por lambda = delta/beta \\
    # 2         -> matriz (-beta/gamma, 0; 0, 1) con referencia los puntos fijos y uno arbitrario
    #
    def simplificar(self):
        paso("Calculamos los puntos fijos")
        fijos = self.puntos_fijos()
        # Identidad: ya simplificada
        if len(fijos) == 0:
            return self
        ref = self._recta.referencia()
        # Elación: matriz de Jordan fija y referencia el punto fijo y dos arbitrarios
        if len(fijos) == 1:
            # Cogemos dos de los puntos de la referencia que no sean el fijo para la nueva referencia
            i = 0 if ref[0] != fijos[0] else 2
            j = 1 if ref[1] != fijos[0] else 2
            lamb = self._matriz[0][1] / -self._matriz[0][0]
            paso("Usamos el punto fijo como primer punto de la nueva referencia y dejamos dos de los otros dos:")
            paso("Cambiamos de referencia theta' = lambda*theta, con lambda = ", \
                self._matriz[0][1], "/", -self._matriz[0][0], "=", lamb)
            # Calculamos la nueva referencia obteniendo
            nueva0 = self._recta[lamb * self._recta.coordenada(fijos[0])]
            nueva1 = self._recta[lamb * self._recta.coordenada(ref[i])]
            nueva2 = self._recta[lamb * self._recta.coordenada(ref[j])]
            paso("R={", nueva0, ",", nueva1, ";", nueva2, "}")
            return homografia_recta(matrix([[1, 1], [0, 1]]), recta_proyectiva(nueva0, nueva1, nueva2))
        # Dos puntos fijos distintos, podría ser involución
        # Aseguramos que el tercer punto de la referencia no coincida con los fijos
        i = 0 if ref[0] != fijos[0] and ref[0] != fijos[1] else \
            1 if ref[1] != fijos[0] and ref[1] != fijos[1] else \
            2
        _no_pasos()
        modulo = self.modulo()
        # Se debe cumplir que {M, N; P, P'} = modulo (con M, N fijos, P arbitrario)
        # Si no, hay que escogerlos al revés
        j = 0 if razon_doble(fijos[0], fijos[1], ref[i], self(ref[i])) == modulo else 1
        _no_pasos(False)
        paso("Usamos los puntos fijos como primeros puntos (en el orden correcto) de la nueva referencia y mantenemos el otro:")
        paso("R={", fijos[j], ",", fijos[(j + 1) % 2], ";", ref[i], "}")
        return homografia_recta(matrix([[modulo, 0], [0, 1]]), recta_proyectiva(fijos[j], fijos[(j + 1) % 2], ref[i]))

    #\m
    # Devuelve esta misma involución en una referencia adecuada de forma que su ecuación quede theta*theta' = 1.
    #
    # Implementación \\
    # Asumiendo que esta homografía es una involución, los dos puntos fijos separan armónicamente cualquier par de la involución.
    # Por tanto, basta con elegir cualquier punto y su imagen como Infinity y 0 y un punto fijo como 1. De esta forma, el otro
    # punto fijo tendrá coordenada -1 y con ello obtenemos la ecuación deseada.
    #
    def simplificar_involucion(self):
        assert self.es_involucion(), "Para simplificar como involucion la homografia debe ser una involucion"
        paso("Calculamos los puntos fijos")
        fijos = self.puntos_fijos()
        ref = self._recta.referencia()
        i = 0 if ref[0] != fijos[0] and ref[0] != fijos[1] else \
            1 if ref[1] != fijos[0] and ref[1] != fijos[1] else \
            2
        paso("Usamos un punto arbitrario y su imagen como primeros puntos de la referencia y uno de los fijos como unidad")
        paso("R={", ref[i], ",", self(ref[i]), ",", fijos[0], "}")
        return homografia_recta(matrix([[0, 1], [1, 0]]), recta_proyectiva(ref[i], self(ref[i]), fijos[0]))

    #\m
    # Devuelve el haz de ecuaciones cuadráticas generado por esta involución.
    #
    # Funciona tanto para puntos (vectores) como parámetros inhomogéneos de la recta.
    #
    # Implementación \\
    # Asumiendo que esta homografía es una involución, calcula las imágenes de los puntos dados. Después, simplemente será
    # necesario componerlos como dos ecuaciones cuyas soluciones sean cada par (theta_i, h(theta_i)).
    #
    # Parámetros \\
    # a: complejo/Infinity/vector(n) - punto del primer par
    # b: complejo/Infinity/vector(n) - punto del segundo par
    #
    def haz_ecuaciones_cuadraticas(self, a, b):
        assert self.es_involucion(), "Para generar el haz de ecuaciones cuadraticas, esta homografia debe ser una involucion"
        assert a != b, "Para generar el haz de ecuaciones cuadraticas se deben especificar puntos distintos"
        _no_pasos()
        if not es_parametro(a):
            a = self._recta.coordenada(a)
        if not es_parametro(b):
            a = self._recta.coordenada(b)
        ap = self(a)
        bp = self(b)
        _no_pasos(False)
        assert ap != b and bp != a, "Para generar el haz de ecuaciones cuadraticas no se puede especificar un par de la involucion"
        paso("Los dos pares especificados son:")
        paso((a, ap), ", ", (b, bp))
        var('mu theta')
        # Cuidado con infinitos
        p0 = (a - theta) if not es_infinito(a) else 1
        p1 = (ap - theta) if not es_infinito(ap) else 1
        p2 = (b - theta) if not es_infinito(b) else 1
        p3 = (bp - theta) if not es_infinito(bp) else 1
        return p0 * p1 + mu * p2 * p3

    #\m
    # Operador *. Devuelve la composición de las homografías ((self o otra)(theta) = self(otra(theta))).
    #
    # Uso: h * k (h y k son homografías de una misma recta).
    #
    # Implementación \\
    # Devuelve una nueva homografía cuya matriz sea el producto de las de ambas.
    #
    # Parámetros \\
    # otra: homografia_recta - homografía con la que componer
    #
    def __mul__(self, otra):
        assert self._recta == otra._recta, "Las homografias deben ser de la misma recta"
        return homografia_recta(self._matriz * otra._matriz, self._recta)

    #\m
    # Operador ^ (ó **). Devuelve el resultado de componer una homografía con sí misma n veces (^-1 devuelve la inversa).
    #
    # Uso h^n (ó h**n) (h es una homografía de la recta y n un entero).
    #
    # Implementación \\
    # Devuelve una nueva homografía cuya matriz es la de esta elevada a n.
    #
    # Parámetros \\
    # n: entero - exponente al que elevar
    #
    def __pow__(self, n):
        return homografia_recta(self._matriz^n, self._recta)

    # Métodos auxiliares

    def __sustituir(self, x):
        # Devuelve el límite para los casos de infinito o división por 0
        return limit(self._expresion_mobius.rhs(), theta = x)

    def __repr__(self):
        return "<Homografia " + str(self._expresion) + " de " + str(self._recta) + ">"

#\c
# Clase que representa una homografía entre dos rectas distintas. En general, se creará utilizando funciones auxiliares.
#
# Se representa mediante una falsa homografía de una recta en sí misma, que para lo único que servirá será para aprovechar sus métodos.
#
class homografia_dos_rectas:

    def __init__(self, homografia, recta1, recta2):
        assert recta1.dimension_ambiente() == recta2.dimension_ambiente(), "Las rectas deben pertenecer al mismo espacio"
        _no_pasos()
        interseccion = recta1.subespacio().interseccion(recta2.subespacio())
        _no_pasos(False)
        assert interseccion.es_punto(), "La interseccion de las rectas debe ser un unico punto"
        self._homografia = homografia
        self._recta_origen = recta1
        self._recta_destino = recta2
        self._interseccion = interseccion.representantes()[0]

    # Métodos accedentes

    #\m
    # Devuelve esta homografía como una aplicacion_proyectiva. Nótese que así dejará de estar restringida a las rectas origen y destino.
    def aplicacion(self):
        return self._homografia.aplicacion()

    #\m
    # Devuelve la matriz asociada a esta homografía.
    def matriz_asociada(self):
        return self._homografia.matriz_asociada()

    #\m
    # Devuelve el punto de interseccion de las rectas de esta homografía.
    def interseccion(self):
        return self._interseccion

    #\m
    # Devuelve la recta de partida de esta homografía.
    def recta_origen(self):
        return self._recta_origen

    #\m
    # Devuelve la referencia en la que está expresada la recta origen.
    def referencia_origen(self):
        return self._recta_origen.referencia()

    #\m
    # Devuelve la recta de destino de esta homografía.
    def recta_destino(self):
        return self._recta_destino

    #\m
    # Devuelve la referencia en la que está expresada la recta destino.
    def referencia_destino(self):
        return self._recta_destino.referencia()

    #\m
    # Devuelve la homografía expresada como una transformación de Möbius.
    def expresion_mobius(self):
        return self._homografia.expresion_mobius()

    #\m
    # Devuelve la ecuación en theta y theta' de esta homografía.
    def ecuacion(self):
        return self._homografia.ecuacion()


    # Otros métodos

    #\m
    # Calcula la imagen mediante esta homografía del punto dado.
    #
    # Uso: h(x) ó h(p) (donde r es una homografia_recta, x un complejo/Infinity y p un punto).
    #
    # Implementación \\
    # Calcula el parámetro no homogéneo asociado a la imagen del punto dado mediante la falsa homografía. \\
    # Si se ha dado un parámetro devolverá el parámetro y si se ha dado un punto devolverá el punto.
    #
    # Parámetros \\
    # x: complejo/Infinity/vector(n) - punto del que se quiere calcular su imagen
    #
    def __call__(self, x):
        # Estamos recibiendo un parámetro no homogéneo, simplemente calculamos su imagen
        if es_parametro(x):
            return self._homografia(x)
        # Asumimos que es un vector, si es del mismo espacio calculamos su coordenada en la recta primero
        assert x in self._recta_origen, "El punto debe pertenecer a la recta origen"
        coord = self._recta_origen.coordenada(x)
        paso("La coordenada de ", x, " en la recta origen es: ", coord, ". Calculamos su imagen y recuperamos el punto")
        return self._recta_destino[self(coord)]

    #\m
    # Calcula el eje de esta homografía utilizando el teorema del eje. Se puede pedir que se devuelva en forma de recta_proyectiva
    # en vez de como subespacio (por defecto). En tal caso, se pueden indicar los puntos que se quiere utilizar para hallar el eje
    # y de los cortes generados se obtendrá la referencia de dicha recta. El primer punto es el generado por ab, el ssegundo ac y el tercero bc.
    # Importante no repetir si se indican los puntos y evitar que uno de ellos sea la intersección (de momento da error).
    #
    # Implementación \\
    # Se escogen los puntos indicados o tres puntos arbitrarios P, Q, R. Se calculan sus imágenes P', Q', R' y se unen las de distinto
    # nombre. La unión de las intersecciones resultará en el eje.
    #
    # Parámetros \\
    # recta: booleano - determina si se quiere obtener un objeto del tipo recta_proyectiva como resultado en vez de un subespacio (False por defecto) \\
    # a: complejo/Infinity - primer punto para obtener el eje (por defecto Infinity) \\
    # b: complejo/Infinity - segundo punto para obtener el eje (por defecto 0) \\
    # c: complejo/Infinity - tercer punto para obtener el eje (por defecto 1)
    #
    def eje(self, recta = False, a = Infinity, b = 0, c = 1):
        assert a != b and a != c and b != c, "Se deben indicar puntos distintos para obtener el eje"
        _no_pasos()
        ap = self(a)
        bp = self(b)
        cp = self(c)
        p1 = self._recta_origen[a]
        p2 = self._recta_destino[ap]
        q1 = self._recta_origen[b]
        q2 = self._recta_destino[bp]
        r1 = self._recta_origen[c]
        r2 = self._recta_destino[cp]
        _no_pasos(False)
        paso("Calculamos las imagenes de tres puntos arbitrarios:", a, ",", b, ",", c)
        paso("Para los puntos indicados se tienen las siguientes transformaciones:")
        paso("P = ", p1, " -> ", p2, " = P'")
        paso("Q = ", q1, " -> ", q2, " = Q'")
        paso("R = ", r1, " -> ", r2, " = R'")
        paso("Unimos los de dsitinto nombre e intersecamos (no se muestran pasos por no llenar la pantalla, pero se pueden calcular facilmente)")
        _no_pasos()
        pq = subespacio(p1, q2)
        pr = subespacio(p1, r2)
        qp = subespacio(q1, p2)
        qr = subespacio(q1, r2)
        rp = subespacio(r1, p2)
        rq = subespacio(r1, q2)
        x1 = pq.interseccion(qp).representantes()[0]
        x2 = pr.interseccion(rp).representantes()[0]
        x3 = qr.interseccion(rq).representantes()[0]
        _no_pasos(False)
        paso("Las intersecciones son:", x1, ", ", x2, ", ", x3)
        if recta:
            return recta_proyectiva(x1, x2, x3)
        return subespacio(x1, x2, x3)

    #\m
    # Determina si esta homografía es una perspectividad.
    #
    # Implementación \\
    # Comprueba que la intersección de las dos rectas sea un punto fijo, esto es, h(O) = a*O, donde h es esta homografía, O el punto
    # de intersección y a un parámetro cualquiera.
    #
    def es_perspectividad(self):
        var('param')
        return len(solve((self(self._interseccion) - param * self._interseccion).list(), param)) == 1

    #\m
    # Calcula el punto centro de esta perspectividad.
    #
    # Implementación \\
    # Calcula el corte de las rectas que forman dos pares arbitrarios que no coincidan con la intersección.
    #
    def centro(self):
        assert self.es_perspectividad(), "El interes de calcular el centro es para perspectividades"
        _no_pasos()
        a = self(Infinity)
        b = self(0)
        c = self(1)
        r1 = subespacio(self._recta_origen[Infinity], self._recta_destino[a])
        r2 = subespacio(self._recta_origen[0], self._recta_destino[b])
        r3 = subespacio(self._recta_origen[1], self._recta_destino[c])
        _no_pasos(False)
        paso("Tomamos los puntos Infinity, 0 y 1 de la recta origen (tres por si alguno fuera la interseccion)")
        paso(Infinity, "->", a)
        paso(0, "->", b)
        paso(1, "->", c)
        paso("Que forman las siguientes rectas y se cortaran en el centro")
        paso(r1)
        paso(r2)
        paso(r3)
        _no_pasos()
        res = r1.interseccion(r2) if not r1.es_punto() and not r2.es_punto() else \
                r1.interseccion(r3) if not r1.es_punto() and not r3.es_punto() else \
                r2.interseccion(r3)
        _no_pasos(False)
        return res.representantes()[0]

    #\m
    # Devuelve la inversa de esta homografía que va de la recta de destino a la recta origen.
    #
    # Implementación \\
    # Devuelve una nueva homografia_dos_rectas cuya falsa homografía sea la inversa de la actual y con las rectas intercambiadas.
    #
    def inversa(self):
        return homografia_dos_rectas(self._homografia^-1, self._recta_destino, self._recta_origen)

    #\m
    # Operador *. Devuelve la composición de las homografías ((self o otra)(theta) = self(otra(theta))). Se debe cumplir que las
    # tres rectas implicadas se corten en un mismo punto o que sólo hay dos.
    #
    # Uso: h * k (h y k son homografías tal que el destino de k es el origen de h en la misma referencia).
    #
    # Implementación \\
    # Devuelve la composición de la falsa homografía de la otra con esta. \\
    # NOTA. Por simplicidad NO se comprueba que las referencias de las rectas que tienen que ser iguales coincidan (aunque debería).
    # Se debe asegurar de todas formas que la referencia de ambas rectas sea la misma o si no el resultado no será válido.
    #
    # Parámetros \\
    # otra: homografia_dos_rectas - homografía con la que componer
    #
    def __mul__(self, otra):
        assert self._recta_origen == otra._recta_destino, \
            "La recta destino de la primera homografia que actua debe ser la misma que el origen de la segunda"
        # Estamos volviendo a la misma recta
        if self._recta_destino == otra._recta_origen:
            return homografia_recta(self._homografia.matriz_asociada() * otra._homografia.matriz_asociada(), self._recta_destino)
        return homografia_dos_rectas(self._homografia * otra._homografia, otra._recta_origen, self._recta_destino)

    #\m
    # Operador ^ (ó **). Sólo aplicable para n = -1. Usar inversa() mejor. Devuelve la inversa de esta homografía.
    #
    # Uso h^-1 (ó h**-1) (h es una homografía entre dos rectas y n un entero).
    #
    # Implementación \\
    # Ver inversa.
    #
    # Parámetros \\
    # n: -1 - exponente al que elevar
    #
    def __pow__(self, n):
        assert n == -1, "No se puede calcular otra cosa que la inversa, puesto que origen y destino no pueden coincidir"
        return self.inversa()

    def __repr__(self):
        return "<Homografia " + str(self._homografia.ecuacion()) + " entre las rectas " + str(self._recta_origen) + " y " + \
            str(self._recta_destino) + ">"

#\f
# Función auxiliar que determina si el parámetro recibido es un punto (tipo vector) o un parámetro
# inhomogéneo de una recta proyectiva.
#
# Parámetros \\
# x: complejo/Infinity/variable/vector(n>1) - elemento a considerar
#
def es_parametro(x):
    # Infinity no tiene método .list(), se ve aparte
    # La lista tendrá un solo elemento a menos que sea una variable (que la trata como un polinomio)
    return es_infinito(x) or len(x.list()) == 1 or type(x) == sage.symbolic.expression.Expression

#\f
# Función auxiliar que determina si el parámetro recibido es infinito o no. Sólo devuelve True para +/-Infinity.
#
# Parámetros \\
# x: complejo/Infinity/variable/vector(n) - elemento a considerar
#
def es_infinito(x):
    return UnsignedInfinityRing(x) is not UnsignedInfinityRing.less_than_infinity()

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
# NOTA. No puede haber cuatro de ellos alineados, pero la función no lo comprueba.
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
    paso("(", ab.implicitas()[0].lhs(), ")(", cd.implicitas()[0].lhs(), ") + lambda * (", ac.implicitas()[0].lhs(), ")(", bd.implicitas()[0].lhs(), ") = 0")
    paso("Por ultimo, forzamos que pase por E =", e, ", y despejamos lambda")
    _no_pasos()
    res = haz.forzar_punto(e)
    _no_pasos(False)
    return res

#\f
# Calcula la razón doble de 4 puntos sobre una cónica, dado un punto base también perteneciente a la misma. \\
# Como una cónica queda definida por 5 puntos no es necesario especificar explícitamente la cónica. \\
# Nótese que la razón doble debería ser la misma sea cual sea el punto base, siempre que se encuentre en la misma cónica
# y no sea uno de los otros.
#
# Ver también métodos similares en la clase parametrizacion_conica.
#
# Implementación \\
# Calcula la razón doble de las rectas que unen la base con cada uno de los puntos (ver razon_doble (recta_proyectiva.sage)).
# Los puntos deben ser todos distintos, aunque no se comprueba.
#
# Parámetros \\
# base: vector(3) - punto base del haz de rectas con el que se va a calcular la razón doble (debe pertenecer a la cónica) \\
# p0: vector(3) - primer punto de la razón doble (debe pertenecer a la cónica) \\
# p1: vector(3) - segundo punto de la razón doble (debe pertenecer a la cónica) \\
# p2: vector(3) - tercer punto de la razón doble (debe pertenecer a la cónica) \\
# p3: vector(3) - cuarto punto de la razón doble (debe pertenecer a la cónica)
#
def razon_doble_conica(base, p0, p1, p2, p3):
    _no_pasos()
    r0 = subespacio(base, p0).dual().punto()
    r1 = subespacio(base, p1).dual().punto()
    r2 = subespacio(base, p2).dual().punto()
    r3 = subespacio(base, p3).dual().punto()
    _no_pasos(False)
    paso("Calculamos las rectas que unen cada uno de los puntos con la base (no se muestra el procedimiento)")
    paso("Los coeficientes de cada una son: ", r0, ", ", r1, ", ", r2, ", ", r3)
    paso("Calculamos su razon doble")
    return razon_doble_puntos(r0, r1, r2, r3)

#\f
# Crea un objeto del tipo homografia_conica dados la parametrización de la cónica y tres puntos y sus imágenes.
# Es válido tanto para puntos como parámetros inhomogéneos. \\
# Si se quiere crear una involución dados dos puntos y sus imágenes se pueden omitir los terceros puntos. La función se podría
# llamar de la siguiente forma: h = crear_homografia_conica(a, ap, b, bp, param = p, involucion = True).
#
# Implementación \\
# Si se dan puntos se obtiene primero su coordenada de la cónica. Una vez obtenidas las coordenadas se crea una nueva
# homografía de una falsa recta en sí misma (ver crear_homografia_recta_theta() en recta_proyectiva.sage), que se utilizará
# para calcular la imagen de un parámetro no homogéneo.
#
# Parámetros \\
# a: complejo/Infinity/vector(3) - primer punto a transformar \\
# ap: complejo/Infinity/vector(3) - imagen del primer punto \\
# b: complejo/Infinity/vector(3) - segundo punto a transformar \\
# bp: complejo/Infinity/vector(3) - imagen del segundo punto \\
# c: complejo/Infinity/vector(3) - tercer punto a transformar (no necesario si se va a crear una involución) \\
# cp: complejo/Infinity/vector(3) - imagen del tercer punto (no necesario si se va a crear una involución) \\
# param: parametrizacion_conica - cónica sobre la que actúa la homografía (por defecto (1, theta, theta^2)) \\
# involucion: booleano - True si se quiere crear una involución, False si no necesariamente (False por defecto)
#
def crear_homografia_conica(a, ap, b, bp, c = None, cp = None, param = None, involucion = False):
    assert c is not None and cp is not None or involucion, "Si no se crea una involucion se deben indicar los terceros puntos"
    if param is None:
        param = parametrizacion_conica(identity_matrix(3))
    if es_parametro(a):
        assert es_parametro(ap) and es_parametro(b) and es_parametro(bp), "Todos los puntos deben ser del mismo tipo"
        # Si no es una involución también hay que comprobar el tercer punto
        if not involucion:
            assert es_parametro(c) and es_parametro(cp), "Todos los puntos deben ser del mismo tipo"
        h = crear_homografia_recta_theta(a, ap, b, bp, c, cp, None, involucion)
        return homografia_conica(h, param)
    theta0 = param.coordenada(a)
    theta0p = param.coordenada(ap)
    theta1 = param.coordenada(b)
    theta1p = param.coordenada(bp)
    theta2 = param.coordenada(c) if not involucion else None
    theta2p = param.coordenada(cp) if not involucion else None
    paso("Las coordenadas en la recta de los puntos que se van a transformar son:")
    paso(theta0, " -> ", theta0p)
    paso(theta1, " -> ", theta1p)
    if not involucion:
        paso(theta2, " -> ", theta2p)
    h = crear_homografia_recta_theta(theta0, theta0p, theta1, theta1p, theta2, theta2p, None, involucion)
    return homografia_conica(h, param)

#\f
# Crea una homografía de una cónica dados un punto y su imagen y el eje de la homografía (como subespacio).
# Es válido tanto para puntos como parámetros inhomogéneos.
#
# Implementación \\
# Obtiene los puntos P, Q de corte del eje con la cónica. Estos serán puntos fijos. Así que se crea una homografía que lleve
# P a P, Q a Q y A a A' (el punto dado). En caso de obtenerse un único punto fijo (el eje era tangente) se calcula la imagen
# del punto extra B dado hallanda el corte X entre BA' con el eje y de nuevo el corte de la recta AX con la cónica, que será
# su imagen. Entonces ya tenemos una homografía que transforma P en P, A en A' y B en B'.
#
# NOTA. Por simplicidad no se hacen comprobaciones sobre puntos coincidentes, pero se debe asegurar que ninguno de los puntos
# que se da como parámetros (tanto a, ap, como extra) coincida con la intersección (o entre ellos). En caso contrario es probable
# encontrarse con resultados extraños o incluso errores.
#
# Parámetros \\
# a: complejo/Infinity/vector(3) - punto a transformar \\
# ap: complejo/Infinity/vector(3) - imagen del punto \\
# param: parametrizacion_conica - cónica sobre la que actúa la homografía \\
# eje: subespacio - eje de la homografía (nótese que se pide un subespacio) \\
# extra: complejo/Infinity/vector(3) - punto extra de la cónica que se utilizará en caso de que el corte del eje con la cónica
# sea un único punto (0)
#
def crear_homografia_conica_eje(a, ap, param, eje, extra = 0):
    if es_parametro(a):
        assert es_parametro(ap), "Todos los puntos deben ser del mismo tipo"
        _no_pasos()
        a = param[a]
        ap = param[ap]
        _no_pasos(False)
    (x, y) = param.interseccion_recta(eje)
    # Distinto punto de corte
    if x != y:
        paso("Los puntos de corte del eje con la conica son: ", x, ", ", y, "; son los puntos fijos")
        paso("Creamos la homografia que ademas transforma: ", a, " -> ", ap)
        return crear_homografia_conica(a, ap, x, x, y, y, param)
    # Mismo punto de corte
    else:
        _no_pasos()
        c = param[extra] if es_parametro(extra) else extra
        _no_pasos(False)
        paso("El punto de corte del eje con las rectas ha sido unico:", x)
        paso("Hallamos la imagen del punto extra escogido:", c, ", uniendolo con A':", ap, " y cortando esta recta con el eje")
        _no_pasos()
        r = subespacio(c, ap)
        z = r.interseccion(eje).punto()
        s = subespacio(a, z)
        (c1, c2) = param.interseccion_recta(s)
        # Elegimos el punto que no sea A
        cp = c2 if matrix([a, c]).rank() == 1 else c1
        _no_pasos(False)
        paso("Resumiendo, la recta corta en el eje en: ", z, ", y uniendolo con A corta la conica en: ", cp)
        return crear_homografia_conica(a, ap, x, x, c, cp, param)

#\f
# Crea una homografía de una cónica dados un punto y su imagen y los puntos fijos. Es válido tanto para puntos como parámetros inhomogéneos.
#
# Implementación \\
# Utiliza la función crear_homografia_conica_eje calculando primero el eje, que será la recta que une los puntos fijos o la
# tangente en el caso de un único punto fijo. Si se tienen dos puntos fijos se puede usar directamente la función principal, sin necesidad
# de calcular el eje.
#
# NOTA. Por simplicidad no se hacen comprobaciones sobre puntos coincidentes, pero se debe asegurar que ninguno de los puntos
# que se da como parámetros (tanto a, ap, como extra) coincida con la intersección (o entre ellos). En caso contrario es probable
# encontrarse con resultados extraños o incluso errores.
#
# Parámetros \\
# a: complejo/Infinity/vector(3) - punto a transformar \\
# ap: complejo/Infinity/vector(3) - imagen del punto \\
# param: parametrizacion_conica - cónica sobre la que actúa la homografía \\
# m: complejo/Infinity/vector(3) - primer punto fijo \\
# n: complejo/Infinity/vector(3) - segundo punto fijo (por defecto ninguno) \\
# extra: complejo/Infinity/vector(3) - punto extra de la cónica que se utilizará si sólo hay un punto fijo (por defecto 0)
#
def crear_homografia_conica_fijos(a, ap, param, m, n = None, extra = 0):
    if es_parametro(m):
        assert n is None or es_parametro(n), "Todos los puntos deben ser del mismo tipo"
        _no_pasos()
        m = param[m]
        if n is not None:
            n = param[n]
        _no_pasos(False)
    paso("El eje une los puntos fijos o es la tangente a la conica")
    eje = subespacio(m, n) if n is not None else param.conica().polar(m)
    paso("El eje tiene ecuacion: ", eje.implicitas()[0])
    return crear_homografia_conica_eje(a, ap, param, eje, extra)

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
        variables = (var('x_var', latex_name = 'x'), var('y_var', latex_name = 'y'), var('z_var', latex_name = 'z'))
        self._vars = vector(variables)
        self._matriz = matriz
        self._ecuacion = (self._vars * matriz * self._vars).factor()

    # Métodos accedentes

    #\m
    # Devuelve la matriz asociada a esta cónica.
    def matriz_asociada(self):
        return self._matriz

    #\m
    # Devuelve la ecuación de esta cónica.
    #
    # Parámetros \\
    # factor: booleano - indica si se quiere la ecuación factorizada (en caso de degeneradas quizá pueda separarla en dos rectas) o no (False por defecto)
    #
    def ecuacion(self, factor = False):
        # Intenta factorizar de nuevo por si acaso
        return self._ecuacion.factor() == 0 if factor else self._ecuacion.expand() == 0

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
    # Devuelve una tupla conteniendo los dos puntos de la cónica dadas dos de sus coordenadas. No se deben dar ni más ni menos.
    # Devolverá dos veces el mismo punto si sólo hay un punto que cumpla esas condiciones.
    #
    # Uso: (P, Q) = c.puntos(x0 = x, y0 = y) (donde c es una cónica, x, y son números); similarmente:
    # P = c.puntos(x0 = x, z0 = z), P = c.puntos(y0 = y, z0 = z)
    #
    # Implementación \\
    # Resuelve la ecuación de la cónica sustituyendo los valores dados para la variable no dada.
    #
    # Parámetros \\
    # x0: complejo - coordenada x del punto que se quiere obtener (si no se especifica será la que se calcule) \\
    # y0: complejo - coordenada y del punto que se quiere obtener (si no se especifica será la que se calcule) \\
    # z0: complejo - coordenada z del punto que se quiere obtener (si no se especifica será la que se calcule)
    #
    def puntos(self, x0 = None, y0 = None, z0 = None):
        assert len(filter(lambda t: t is None, [x0, y0, z0])) == 1, "Se debe especificar el valor de dos variables exactamente"
        if x0 is None:
            sol = map(lambda s: s.rhs(), solve(self.ecuacion().substitute(self._vars[1] == y0, self._vars[2] == z0), self._vars[0]))
            return (vector([sol[0], y0, z0]), vector([sol[-1], y0, z0]))
        if y0 is None:
            sol = map(lambda s: s.rhs(), solve(self.ecuacion().substitute(self._vars[0] == x0, self._vars[2] == z0), self._vars[1]))
            return (vector([x0, sol[0], z0]), vector([x0, sol[-1], z0]))
        if z0 is None:
            sol = map(lambda s: s.rhs(), solve(self.ecuacion().substitute(self._vars[0] == x0, self._vars[1] == y0), self._vars[2]))
            return (vector([x0, y0, sol[0]]), vector([x0, y0, sol[-1]]))

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
    # Devuelve una nueva cónica expresada en la nueva referencia R', dada la matriz del cambio
    # de referencia de R a R'.
    #
    # La matriz del cambio de referencia debería proceder de cambiar_referencia (espacios.sage).
    #
    # Implementación \\
    # Si A es la matriz asociada a esta cónica y M la matriz del cambio de referencia entre dos referencias
    # (que esta cónica desconoce), se crea una nueva cónica cuya matriz asociada sea M^-t * A * M^-1.
    #
    # Parámetros \\
    # matriz_cambio: matriz(3, 3) - matriz que representa el cambio de referencia: x' = Ax
    #
    def cambiar_referencia(self, matriz_cambio):
        assert matriz_cambio.nrows() == 3, "La matriz del cambio de referencia de una cónica es 3x3"
        assert matriz_cambio.det() != 0, "Una matriz de cambio de referencia debe ser invertible"
        paso("La matriz de la conica en la nueva referencia se calcula:", (matriz_cambio.T)^-1, self._matriz, matriz_cambio^-1)
        return conica((matriz_cambio.T)^-1 * self._matriz * matriz_cambio^-1)

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
    # Devuelve la cónica degenerada tangente a esta cónica que pasa por P.
    #
    # Implementación \\
    # Devuelve la cónica de matriz A * p * p^t * A - (p^t * A * p) * A, donde A es la matriz asociada a esta cónica.
    #
    # Parámetros \\
    # p: vector(3) - punto del que obtener las tangentes
    #
    def tangentes(self, p):
        x = matrix([p]).T
        paso("La conica degenerada tangente a esta que pasa por el punto es la de matriz")
        paso(self._matriz, x, x.T, self._matriz, " - (", x.T, self._matriz, x, ")", self._matriz)
        return conica(self._matriz * x * x.T * self._matriz - (x.T * self._matriz * x)[0][0] * self._matriz)

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
    # Devuelve una tupla conteniendo los dos puntos de intersección de una recta (expresada como recta_proyectiva)
    # con esta cónica. Si la recta fuera tangente (esto es, corta en un punto doble) devolverá dos veces el mismo punto.
    #
    # Implementación \\
    # Tomado un punto genérico de la recta, P(theta), resuelve P^t(theta) * A * P(theta) = 0, siendo A la matriz que
    # representa esta cónica.
    #
    # Parámetros \\
    # r: recta_proyectiva(dim=2) - recta con la que intersecar
    #
    def interseccion_recta(self, r):
        assert r.dimension_ambiente() == 2, "La recta debe ser del plano"
        var('theta', latex_name = '\\theta')
        p = r[theta]
        ec = (p * self._matriz * p).expand()
        assert ec != 0, "La recta con la que se intenta intersecar esta contenida en la conica"
        # Si queda una ecuación de primer grado o grado cero es que una solución era infinito
        sol = [] if ec.degree(theta) == 2 else [theta == Infinity]
        paso("Resolvemos: ", matrix([p]), self._matriz, matrix([p]).T, " = ", ec, " = 0")
        sol = sol + solve(ec, theta)
        paso(sol)
        paso("Sustituimos en la recta para obtener los dos puntos")
        return (r[sol[0].rhs()], r[sol[-1].rhs()])

    #\m
    # Devuelve una tupla con los cuatro puntos de intersección de esta cónica con otra cónica, dada su parametrización.
    #
    # NOTA. No será posible intersecar cónicas ambas degeneradas, pues ninguna se podría parametrizar. En tal caso,
    # es recomendable buscar las rectas que forman alguna de ellas e intersecar con las rectas.
    #
    # Implementación \\
    # Resuelve la ecuación c^t(theta) * A * c(theta), donde c es la parametrización de la otra cónica y A es la matriz
    # de esta.
    #
    # Parámetros \\
    # param: parametrizacion_conica - parametrización de la cónica con la que intersecar
    #
    def interseccion_conica(self, param):
        theta = var('theta_var', latex_name = '\\theta')
        _no_pasos()
        c = param[theta]
        _no_pasos(False)
        ec = (c * self._matriz * c).expand()
        paso("Resolvemos: ", matrix([c]), self._matriz, matrix([c]).T, " = ", ec, " = 0")
        assert ec != 0, "La conica con la que se intenta intersecar esta contenida en la conica"
        # Por cada término que falte hay una solución que es infinito
        inf = [theta == Infinity for i in range(4, ec.degree(theta), -1)]
        # Resolvemos con multiplicidad
        (res, mult) = solve(ec, theta, multiplicities = True)
        # Ponemos todo junto
        sol = inf + flatten([[res[i]] * mult[i] for i in range(len(res))])
        paso(sol)
        return tuple(map(lambda t: param[t.rhs()], sol))

    #\m
    # Devuelve el tercer punto que forma un triángulo autopolar junto con los otros dos dados, asumiendo que cada uno
    # está en la polar del otro.
    #
    # Si se tiene sólo un punto se debe calcular antes su polar y de ahí elegir un punto arbitrario para que sea
    # el segundo.
    #
    # Implementación \\
    # Interseca las rectas polares de los puntos dados.
    #
    # Parámetros \\
    # p: vector(3) - un punto del triángulo \\
    # q: vector(3) - otro punto del triángulo
    #
    def triangulo_autopolar(self, p, q):
        paso("Polar del primer punto")
        s = self.polar(p)
        assert q in s, "Cada punto debe pertenecer a la polar del otro"
        paso("Polar del segundo punto")
        t = self.polar(q)
        paso("Intersecamos")
        assert p in t, "Cada punto debe pertenecer a la polar del otro"
        return s.interseccion(t).punto()

    #\m
    # Determina si P, Q, R forman un triángulo autopolar.
    #
    # Implementación \\
    # Comprueba que Q y R pertenezcan a la polar de P.
    #
    # Parámetros \\
    # p: vector(3) - primer punto del triángulo autopolar \\
    # q: vector(3) - segundo punto del triángulo autopolar \\
    # r: vector(3) - tercer punto del triángulo autopolar
    #
    def es_triangulo_autopolar(self, p, q, r):
        _no_pasos()
        polar = self.polar(p)
        _no_pasos(False)
        return q in polar and r in polar

    #\m
    # Devuelve el punto unidad E para la referencia {P, Q, R; E} que hace que esta cónica tenga ecuación x^2 + y^2 + z^2 = 0
    # en complejos, quizá algunos negativos en reales. Nótese que no cambia la cónica, sino que se debe obtener primero la
    # matriz del cambio de referencia (ver cambiar_referencia() en espacios.sage) y luego cambiar la referencia de esta
    # cónica.
    #
    # NOTA. Los puntos P, Q, R deben formar un triángulo autopolar.
    #
    # Implementación \\
    # Con un punto E arbitrario la ecuación de esta cónica es ax^2 + by^2 + cz^2 = 0, donde a, b y c se calculan mediante
    # B(P, P), B(Q, Q), B(R, R), respectivamente (donde B es la forma bilineal asociada a esta cónica). Entonces sólo necesitamos
    # eliminar estos a, b y c. Esto se hace eligiendo los representantes de P, Q y R de forma que E = [beta1 * P + beta2 * Q + beta3 * R],
    # donde a*beta1^2 = 1, b*beta2^2 = 1 y c*beta3^2 = 1 (si a, b o c son 0 se elige betai = 1). \\
    # Para el caso real se igualan a -1 si no se puede resolver la ecuación.
    #
    # Parámetros \\
    # p: vector(3) - primer punto del triángulo autopolar \\
    # q: vector(3) - segundo punto del triángulo autopolar \\
    # r: vector(3) - tercer punto del triángulo autopolar \\
    # real: booleano - determina si la cónica debe ser necesariamente real o no (False por defecto)
    #
    def referencia_autopolar(self, p, q, r, real = False):
        assert self.es_triangulo_autopolar(p, q, r), "Las referencias autopolares se forman con triangulos autopolares"
        a = self(p)
        b = self(q)
        c = self(r)
        paso("a = B(p,p) = ", a, "; b = B(q,q) = ", b, "; c = B(r,r) = ", c)
        paso("Calculamos betai para que a*beta1 = a*beta1^2 = 1, b*beta2^2 = 1 y c*beta3^2 = 1 (o betai = 1 si su respectivo es 0)")
        if not real:
            beta1 = sqrt(1/a).simplify_full() if a != 0 else 1
            beta2 = sqrt(1/b).simplify_full() if b != 0 else 1
            beta3 = sqrt(1/c).simplify_full() if c != 0 else 1
            paso("beta1 = ", beta1, "; beta2 = ", beta2, "; beta3 = ", beta3)
            return beta1 * p + beta2 * q + beta3 * r
        paso("Si alguna ecuacion no se puede resolver en reales se iguala a -1")
        beta1 = sqrt(1/abs(a)).simplify_full() if a != 0 else 1
        beta2 = sqrt(1/abs(b)).simplify_full() if b != 0 else 1
        beta3 = sqrt(1/abs(c)).simplify_full() if c != 0 else 1
        paso("beta1 = ", beta1, "; beta2 = ", beta2, "; beta3 = ", beta3)
        return beta1 * p + beta2 * q + beta3 * r

    #\m
    # Devuelve una parametrización de la cónica, dados tres puntos contenidos en ella que actuarán como primero, tercero y unidad de
    # la referencia en la que la cónica tiene ecuación y^2=xz.
    #
    # Implementación \\
    # Calcula las polares (tangentes por estar en la cónica) de los dos primeros puntos: P, R, para obtener el tercero, que será el corte
    # de ambas rectas: Q. Tomando como referencia entonces {P, Q, R; E}, la ecuación de la cónica sería y^2=xz. Por tanto, se puede
    # parametrizar como (1, theta, theta^2).
    #
    # Parámetros \\
    # p: vector(3) - primer punto de la referencia (debe pertenecer a la cóncia) \\
    # r: vector(3) - tercer punto de la referencia (debe pertenecer a la cóncia) \\
    # e: vector(3) - punto unidad de la referencia (debe pertenecer a la cóncia)
    #
    def parametrizar(self, p, r, e):
        assert not self.es_degenerada(), "No se pueden parametrizar conicas degeneradas"
        assert p in self and r in self and e in self, "Los puntos deben pertenecer a la conica"
        paso("Calculamos las polares de: ", p, " y ", r)
        polarp = self.polar(p)
        polarr = self.polar(r)
        _no_pasos()
        q = polarp.interseccion(polarr).punto()
        _no_pasos(False)
        paso("Intersecamos: ", q)
        paso("Calculamos la matriz asociada a la referencia  {", p, ", ", q, ", ", r, "; ", e, "} para los cambios de referencia futuros")
        m = matriz_asociada(matrix([p, q, r, e]).T)
        return parametrizacion_conica(identity_matrix(3), m)

    # ELIMINADO DE MOMENTO PORQUE PARECE NO FUNCIONAR BIEN
    '''
    #
    # Devuelve una tupla conteniendo las dos rectas (como subespacios) en que se descompone esta cónica degenerada.
    # Si la cónica es una recta doble ambas rectas serán la misma.
    #
    # NOTA. Aquí se utiliza el procedimiento que se describe a continuación porque es más sencillo de programar, pero hay otros métodos, como el de resolver
    # la ecuación en alguna de las variables, o hallar el punto de intersección de las rectas y luego cortar la cónica con rectas sencillas, como x=0, y=0 o z=0
    # para hallar los otros dos puntos y luego unirlos con el primero. Además, simplemente accediendo a la ecuación de la cónica con factor = True (ver ecuacion())
    # se puede ver la cónica factorizada, pero sin ningún tipo de procedimiento especial. Esto está por si fuera necesario.
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
    '''

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
        paso(m)
        paso(uxv1)
        paso(uxv2)
        paso(uxv3)
        m01 = (uxv1 * uxv2).simplify_full()
        m02 = (uxv1 * uxv3).simplify_full()
        m12 = (uxv2 * uxv3).simplify_full()
        s1 = 1 if m01 == m[0][1] else -1
        s2 = 1 if m02 == m[0][2] else -1
        s3 = 1 if m12 == m[1][2] else -1
        return (s1 * s2 * uxv1, s1 * s3 * uxv2, s2 * s3 * uxv3)

    def _tipo(self):
        return 0

    def __repr__(self):
        return "<Conica de ecuacion " + str(self.ecuacion()) + ">"

#\c
# Clase que representa una parametrización de una cónica no degenerada de la forma (p1(theta), p2(theta), p3(theta)),
# siendo pi polinomios de grado 2.
#
class parametrizacion_conica:

    #\i
    # Construye una parametrización de una cónica dadas la matriz M que determina las coordenadas de cada punto, que se
    # obtienen como M * (1, theta, theta^2). M será la matriz que hace que la cónica se exprese como y^2=xz.
    #
    # También se puede indicar una matriz de cambio de referencia extra para casos del tipo (por ejemplo) parametrizar una
    # cónica como (1, theta, theta^2).
    #
    # Parámetros \\
    # matriz: matriz(3, 3) - matriz que representa los puntos genéricos de la parametrización como se explica arriba \\
    # matriz_cambio: matriz(3, 3) - matriz del cambio de referencia para convertir la referencia en que se encuentra la
    # parametrización a la canónica (identidad por defecto, esto es, ya se expresa en la canónica)
    #
    def __init__(self, matriz, matriz_cambio = identity_matrix(3)):
        assert matriz.dimensions() == (3, 3) and matriz_cambio.dimensions() == (3, 3), \
            "Las matrices de la parametrizacion y cambio de referencia deben ser 3x3"
        assert matriz.rank() == 3 and matriz_cambio.rank() == 3, \
            "Las matrices de la parametrizacion y cambio de referencia deben ser invertibles"
        self._matriz = matriz
        self._matriz_cambio = matriz_cambio

    # Métodos accedentes

    #\m
    # Devuelve un punto genérico de la parametrización expresado en la referencia en que se encuentra la parametrización.
    #
    # Diferenciar del operador [] en que c[theta_0] SÍ que cambia la referencia a la referencia original. Este método NO.
    #
    def punto_generico(self):
        theta = var('theta_var', latex_name = '\\theta')
        return self._matriz * vector([1, theta, theta^2])

    #\m
    # Devuelve el punto con coordenada theta.
    #
    # Uso: c[theta_0] (theta_0 es la coordenada del punto y c es parametrizacion_conica).
    #
    # Si se quiere un punto genérico de la cónica, sustituir por una variable (ej: c[x] devuelve un punto dependiente de x).
    #
    # Parámetros \\
    # theta: complejo/Infinity - coordenada del punto de la cónica que se quiere obtener
    #
    def __getitem__(self, theta):
        y = self.__punto_ref(theta)
        if self._matriz_cambio != identity_matrix(3):
            paso("En la referencia de la parametrizacion el punto de coordenada: ", theta, " es: ", y)
            paso("Cambiamos de referencia con la matriz del cambio: ", self._matriz_cambio)
        return self._matriz_cambio * y

    #\m
    # Devuelve la referencia en que está expresada esta parametrización en función de la referencia canónica
    # (es decir, las columnas de la matriz del cambio).
    def referencia(self):
        return self._matriz_cambio.columns()

    #\m
    # Devuelve la referencia en que esta cónica tiene coordenadas (1, theta, theta^2) (las columnas de la matriz que
    # representa los puntos genéricos de la parametrización).
    #
    # NOTA. No estoy muy seguro pero creo que esto está dado en función de la primera referencia.
    def referencia2(self):
        return self._matriz.columns()

    # Otros métodos

    #\m
    # Devuelve la cónica que representa esta parametrización.
    #
    # Implementación \\
    # Conociendo la matriz que representa la parametrización y la del cambio de referencia (ambas son cambios de
    # referencia al fin y al cabo), la cónica tendrá ecuación y^2=xz, por tanto, hay que cambiar de referencia
    # esa cónica primero por la matriz que representa la parametrización, P, y luego por la del cambio de referencia,
    # M: M^-t * P^-t * A * P^-1 * M^-1
    #
    def conica(self):
        a = matrix([[0, 0, -1], [0, 2, 0], [-1, 0, 0]])
        paso("La matriz de la conica que representa esta parametrizacion sera:")
        paso(self._matriz_cambio, "^-t * ", self._matriz, "^-t * ", a, " * ", self._matriz, "^-1 * ", self._matriz_cambio, "^-1")
        return conica(self._matriz_cambio.T^-1 * self._matriz.T^-1 * a * self._matriz^-1 * self._matriz_cambio^-1)

    #\m
    # Devuelve la coordenada no homogénea asociada al punto P dado.
    #
    # Implementación \\
    # Resuelve la ecuación (p1(theta) : p2(theta) : p3(theta)) == P. Se asume que el punto está en la cónica.
    #
    # Parámetros \\
    # p: vector(3) - punto de la cónica del que se quiere conocer su coordenada
    #
    def coordenada(self, p):
        assert len(p) == 3, "El punto debe ser del plano"
        t = var('theta_var', latex_name = '\\theta')
        l = var('lambda_var', latex_name = '\\lambda')
        x = self._matriz * vector([1, t, t^2])
        q = self._matriz_cambio * x
        lp = l * p
        paso("Resolvemos: ", matrix([q]).T, " = ", self._matriz_cambio, matrix([x]).T, " = ", matrix([lp]).T)
        sol = solve((q - lp).list(), t, l)
        paso(sol)
        # No hay solución
        if len(sol) == 0:
            _no_pasos()
            r = matrix([self[Infinity], p]).rank() == 1
            _no_pasos(False)
            assert r, "El punto debe pertenecer a la conica"
            paso("No habia solucion, pero el punto pertenecia a la conica, luego su coordenada es ", Infinity)
            return Infinity
        return sol[0][0].rhs()

    #\m
    # Devuelve una tupla conteniendo los dos puntos de intersección de una recta (expresada como subespacio)
    # con esta cónica. Si la recta fuera tangente (esto es, corta en un punto doble) devolverá dos veces el mismo punto.
    #
    # La recta puede venir dada tanto para coordenadas en la referencia de la recta como para coordenadas en la referencia
    # original. Esto sólo cambiará el orden del producto.
    #
    # Implementación \\
    # Tomado un punto genérico de la cónica, P(theta), lo sustituye en la ecuación implícita de la recta y la resuelve.
    #
    # Parámetros \\
    # r: subespacio(dim=1, dim_ambiente=2) - recta con la que intersecar \\
    # original: booleano - determina si se usa la referencia original o no (True por defecto)
    #
    def interseccion_recta(self, r, original = True):
        assert r.dim() == 1, "El subespacio debe ser una recta"
        assert r.dimension_ambiente() == 2, "La recta debe ser del plano"
        theta = var('theta_var', latex_name = '\\theta')
        p = self._matriz * vector([1, theta, theta^2])
        if original:
            p = self._matriz_cambio * p
        _no_pasos()
        v = r.dual().punto()
        _no_pasos(False)
        ec = (v * p).expand()
        paso("Resolvemos: ", ec, " = 0")
        assert ec != 0, "La recta con la que se intenta intersecar esta contenida en la conica"
        sol = [] if ec.degree(theta) == 2 else [theta == Infinity]
        sol = sol + solve(ec, theta)
        paso(sol)
        if not original:
            paso("Los puntos obtenidos han sido: ", self.__punto_ref(sol[0].rhs()), " y ", self.__punto_ref(sol[-1].rhs()))
            paso("Cambiamos de referencia con la matriz del cambio: ", self._matriz_cambio)
        res1 = self[sol[0].rhs()] if original else self._matriz_cambio * self.__punto_ref(sol[0].rhs())
        res2 = self[sol[-1].rhs()] if original else self._matriz_cambio * self.__punto_ref(sol[-1].rhs())
        return (res1, res2)

    #\m
    # Operador in. Determina si un punto está contenido en esta cónica o no. El punto dado debe estar en la
    # referencia inicial. Es decir, si la parametrización es (1 : theta : theta^2), por ejemplo, el punto (1:1:1)
    # no tiene por qué estar contenido en la cónica si la referencia es otra.
    #
    # Uso: P in c (P es un punto y c una parametrización de una cónica).
    #
    # Implementación \\
    # Comprueba si la ecuación (p1(theta) : p2(theta) : p3(theta)) == P tiene solución.
    #
    # Parámetros \\
    # punto: vector(3) - punto que comprobar si pertenece a la cónica
    #
    def __contains__(self, punto):
        assert len(punto) == 3, "El punto debe ser del plano"
        t = var('theta_var', latex_name = '\\theta')
        l = var('lambda_var', latex_name = '\\lambda')
        x = self._matriz_cambio * self._matriz * vector([1, t, t^2])
        lp = l * punto
        sol = solve((x - lp).list(), t, l)
        # No hay solución, pero quizá sea el punto infinito
        if len(sol) == 0:
            _no_pasos()
            r = matrix([self[Infinity], punto]).rank() == 1
            _no_pasos(False)
            return r
        return True

    #\f
    # Calcula la razón doble {A, B; C, D}. Es válido tanto para puntos como para parámetros inhomogéneos
    # de la cónica.
    #
    # Implementación \\
    # Una vez obtenida la coordenada de cada punto se utiliza razon_doble_theta (ver recta_proyectiva.sage).
    #
    # Parámetros \\
    # a: complejo/Infinity/vector(3) - primer punto de la razón doble \\
    # b: complejo/Infinity/vector(3) - segundo punto de la razón doble \\
    # c: complejo/Infinity/vector(3) - tercer punto de la razón doble \\
    # d: complejo/Infinity/vector(3) - cuarto punto de la razón doble
    #
    def razon_doble(self, a, b, c, d):
        if es_parametro(a):
            assert es_parametro(b) and es_parametro(c) and es_parametro(d), "Todos los puntos deben ser del mismo tipo"
            return razon_doble_theta(a, b, c, d)
        _no_pasos()
        theta0 = self.coordenada(a)
        theta1 = self.coordenada(b)
        theta2 = self.coordenada(c)
        theta3 = self.coordenada(d)
        _no_pasos(False)
        paso("Las coordenadas de los puntos son: ", theta0, ", ", theta1, ", ", theta2, ", ", theta3, \
             ", respectivamente (no se muestra procedimiento)")
        return razon_doble_theta(theta0, theta1, theta2, theta3)

    #\f
    # Determina si los puntos dados forman una cuaterna armónica, esto es, {A, B; C, D} = -1. Es válido tanto
    # para puntos como para parámetros inhomogéneos de una cónica.
    #
    # Parámetros \\
    # a: complejo/Infinity/vector(3) - primer punto de la razón doble \\
    # b: complejo/Infinity/vector(3) - segundo punto de la razón doble \\
    # c: complejo/Infinity/vector(3) - tercer punto de la razón doble \\
    # d: complejo/Infinity/vector(3) - cuarto punto de la razón doble
    #
    def es_cuaterna_armonica(self, a, b, c, d):
        return bool(self.razon_doble(a, b, c, d) == -1)

    #\f
    # Devuelve el conjugado armónico de c respecto de a y b. Es válido tanto para puntos como para parámetros inhomogéneos
    # de una cónica.
    #
    # Parámetros \\
    # a: complejo/Infinity/vector(3) - primer punto de la razón doble \\
    # b: complejo/Infinity/vector(3) - segundo punto de la razón doble \\
    # c: complejo/Infinity/vector(3) - tercer punto de la razón doble
    #
    def conjugado_armonico(self, a, b, c):
        if es_parametro(a):
            assert es_parametro(b) and es_parametro(c), "Todos los puntos deben ser del mismo tipo"
            return conjugado_armonico_theta(a, b, c)
        _no_pasos()
        theta0 = self.coordenada(a)
        theta1 = self.coordenada(b)
        theta2 = self.coordenada(c)
        _no_pasos(False)
        paso("Las coordenadas de los puntos son: ", theta0, ", ", theta1, ", ", theta2, ", respectivamente (no se muestra procedimiento)")
        return self[conjugado_armonico_theta(theta0, theta1, theta2)]

    # Métodos auxiliares

    def __punto_ref(self, theta):
        x = vector([1, theta, theta^2]) if not es_infinito(theta) else vector([0, 0, 1])
        return self._matriz * x

    def _tipo(self):
        return 1

    def __repr__(self):
        t = var('theta_var', latex_name = '\\theta')
        return "<Conica parametrizada como " + str(self._matriz * vector([1, t, t^2])) + " en la referencia " + str(self.referencia()) + ">"

#\c
# Clase que representa una homografía de una cónica en sí misma.
#
# Se representa mediante una falsa homografía de una recta en sí misma, que para lo único que servirá será para aprovechar sus métodos.
#
class homografia_conica:

    #\i
    # Construye una homografía de una cónica dada una homografía de una falsa recta en sí misma y la parametrización de la cónica.
    #
    # En general, no se usará este constructor directamente, sino que se creará a partir de funciones creadoras externas.
    #
    # Parámetros \\
    # homografia: homografia_recta - homografía que convierte la coordenada inhomogénea de un punto de la cónica en su imagen \\
    # param: parametrizacion_conica - conica sobre la que actúa la homografía
    #
    def __init__(self, homografia, param):
        self._homografia = homografia
        self._conica = param

    # Métodos accedentes

    #\m
    # Devuelve la matriz asociada a esta homografía.
    def matriz_asociada(self):
        return self._homografia.matriz_asociada()

    #\m
    # Devuelve la parametrización de la cónica de esta homografía.
    def parametrizacion_conica(self):
        return self._conica

    #\m
    # Devuelve la referencia en la que está expresada la cónica.
    def referencia(self):
        return self._conica.referencia()

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
    # Uso: h(x) ó h(p) (donde r es una homografia_conica, x un complejo/Infinity y p un punto).
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
        # Asumimos que es un vector
        assert x in self._conica, "El punto debe pertenecer a la cónica"
        _no_pasos()
        coord = self._conica.coordenada(x)
        _no_pasos(False)
        paso("La coordenada de ", x, " en la conica es: ", coord, ". Calculamos su imagen y recuperamos el punto")
        return self._conica[self(coord)]

    #\m
    # Devuelve los parámetros de los puntos fijos de esta homografía. Devuelve una lista vacía en caso de que sea la identidad.
    #
    # Implementación \\
    # Ver homografia_recta.puntos_fijos_theta() (recta_proyectiva.sage).
    #
    def puntos_fijos_theta(self):
        return self._homografia.puntos_fijos_theta()

    #\m
    # Devuelve los puntos fijos de esta homografía. Devuelve una lista vacía en caso de que sea la identidad.
    #
    # Implementación \\
    # Ver homografia_recta.puntos_fijos() (recta_proyectiva.sage).
    #
    def puntos_fijos(self):
        fijos = self.puntos_fijos_theta()
        paso("Las coordenadas de los puntos fijos son: ", fijos)
        paso("Recuperamos los puntos:")
        return map(lambda theta: self._conica[theta], fijos)

    #\m
    # Calcula el eje de esta homografía.
    #
    # Implementación \\
    # Se unen los puntos fijos en caso de haber dos o se obtiene la recta tangente al punto fijo en caso de haber uno.
    #
    def eje(self):
        paso("Calculamos los puntos fijos:")
        fijos = self.puntos_fijos()
        assert len(fijos) != 0, "No tiene sentido calcular el eje de la identidad"
        paso("Los puntos fijos son: ", fijos)
        if len(fijos) == 2:
            paso("Los unimos")
            return subespacio(fijos[0], fijos[1])
        paso("Tomamos la recta tangente (para ello necesitamos saber la matriz de la conica)")
        return self._conica.conica().polar(fijos[0])

    #\m
    # Determina si esta homografía es una involución.
    #
    # Implementación \\
    # Comprueba que la homografía al cuadrado sea la identidad.
    #
    def es_involucion(self):
        h2 = (self._homografia^2).matriz_asociada()
        paso("Comprobamos que esta homografia al cuadrado sea la identidad: h^2 = ", h2)
        return h2 / max(max(h2, key = abs), key = abs) == identity_matrix(2)

    #\m
    # Calcula el módulo de esta homografía. Ver homografia_recta.modulo() (recta_proyectiva.sage).
    #
    def modulo(self):
        return self._homografia.modulo()

    #\m
    # Calcula el punto centro de esta involución.
    #
    # Implementación \\
    # Calcula el corte de las rectas tangentes de los puntos fijos. Se asume que esta homografía es una involución.
    #
    def centro(self):
        _no_pasos()
        inv = self.es_involucion()
        _no_pasos(False)
        assert inv, "El interes de calcular el centro es para involuciones"
        _no_pasos()
        fijos = self.puntos_fijos()
        _no_pasos(False)
        paso("Calculamos los puntos fijos: ", fijos, "; y ahora necesitaremos la ecuacion de la conica")
        c = self._conica.conica()
        paso(c.ecuacion())
        paso("Calculamos las tangentes (o polares porque pertenecen a ella) de los puntos a la conica")
        t1 = c.polar(fijos[0])
        t2 = c.polar(fijos[1])
        paso("Intersecamos")
        return t1.interseccion(t2).punto()

    #\m
    # Descompone esta homografía en dos involuciones dado el centro de la primera (s2), de forma que h = s2*s1, donde h es esta homografía.
    # Téngase en cuenta la forma de calcular las imágenes: primero se aplica s1 y luego s2 (s2(s1(x))).
    #
    # El centro dado debe encontrarse en el eje de la homografía, porque si no se complica demasiado.
    #
    # El resultado se da como una tupla que contiene dos objetos del tipo involucion_conica.
    #
    # Implementación \\
    # Crea una primera involución con el centro dado. Se obtienen sus puntos fijos. Uno de ellos al menos no es fijo en la homografía.
    # Este punto fijo de la involución se une con una recta con su imagen original y se interseca con el eje de la homografía. Este será
    # el centro de la segunda homografía (si tenía un solo punto fijo ambas involuciones lo mantienen y si tenía dos ambas involuciones
    # los intercambian).
    #
    # Parámetros \\
    # centro: vector(3) - punto que actuará como centro de la primera involución (debe pertenecer al eje de la homografía) \\
    # tipo_conica: booleano - determina si la involución se creará con la parametrización de la cónica o la cónica como tal
    # (False por defecto; no se mostrará procedimiento de la transformación igualmente)
    #
    def descomponer(self, centro, tipo_conica = False):
        paso("Necesitaremos ele eje, lo calculamos:")
        eje = self.eje()
        assert centro in eje, "El centro dado para descomponer en involuciones debe pertenecer al eje"
        _no_pasos()
        c = self._conica if not tipo_conica else self._conica.conica()
        _no_pasos(False)
        s1 = involucion_conica(c, centro)
        paso("Obtenemos los puntos fijos de la involucion dada por el centro: ", centro, "; y nos quedamos con uno que no coincida con el de esta")
        (m, n) = s1.puntos_fijos()
        _no_pasos()
        fijos = self.puntos_fijos()
        # Si la homografía tiene dos fijos distintos cogemos cualquiera; si no, el que sea distinto
        x = m if len(fijos) == 2 or matrix([m, fijos[0]]).rank() == 2 else n
        hx = self(x)
        _no_pasos(False)
        paso("Calculamos la imagen a traves de esta homografia del punto fijo de la involucion elegido: h", x, " = ", hx)
        paso("Los unimos e intersecamos con el eje para obtener el otro centro:")
        v = subespacio(x, hx).interseccion(eje).punto()
        s2 = involucion_conica(c, v)
        return (s1, s2)

    #\m
    # Operador *. Devuelve la composición de las homografías ((self o otra)(theta) = self(otra(theta))). Ambas deben ser de la misma
    # cónica, con la misma parametrización.
    #
    # Uso: h * k (h y k son homografías de la misma cónica).
    #
    # Implementación \\
    # Crea una nueva homografía de la cónica cuya homografía de la recta asociada sea el producto de la homografía asociada a esta cónica
    # por la de la otra. Ver operador * (__mul__) de homografia_recta (recta_proyectiva.sage). \\
    # NOTA. Por simplicidad NO se comprueba que las parametrizaciones sean iguales. En caso contrario, el resultado no será válido
    # (nótese que se usa la cónica de esta homografía para crear la nueva).
    #
    # Parámetros \\
    # otra: homografia_conica - homografía con la que componer
    #
    def __mul__(self, otra):
        return homografia_conica(self._homografia * otra._homografia, self._conica)

    #\m
    # Operador ^ (ó **). Devuelve el resultado de componer una homografía con sí misma n veces (^-1 devuelve la inversa).
    #
    # Uso h^n (ó h**n) (h es una homografía de una cónica y n un entero).
    #
    # Implementación \\
    # Crea una nueva homografía de la cónica cuya homografía de la recta asociada sea la homografía asociada a esta cónica elevada a n.
    # Ver operador ^/** (__pow__) de homografia_recta (recta_proyectiva.sage).
    #
    # Parámetros \\
    # n: entero - exponente al que elevar
    #
    def __pow__(self, n):
        return homografia_conica(self._homografia^n, self._conica)

    def __repr__(self):
        return "<Homografia " + str(self._homografia.ecuacion()) + "de la conica " + str(self._conica) + ">"

#\c
# Clase que define una involución sobre una cónica, conocido alguno de sus elementos.
#
class involucion_conica:

    #\i
    # Construye una involución sabiendo la cónica al menos unos de sus elementos: vértice, puntos fijos o eje.
    #
    # Nótese que si se da una cónica y el eje, este debe ser del tipo recta_proyectiva y si se da una parametrización de la cónica
    # el eje debe ser un subespacio.
    #
    # En caso de dar más información de la necesaria, sólo se hará caso al parámetro que aparezca antes. \\
    # Por ejemplo involucion_conica(c, v, e) ignorará el eje y sólo se quedará con el vértice; o
    # involucion_conica(c, eje = e, fijos = (m, n)) ignorará los puntos fijos.
    #
    # Implementación \\
    # Si no se da el vértice se calcula como el polo del eje: la recta que une los puntos fijos.
    #
    # Parámetros \\
    # conica: conica/parametrizacion_conica - cónica sobre la que actúa la involución \\
    # vertice: vector(3) - vértice de la involución (por defecto el cálculo que se explica más arriba) \\
    # eje: recta_proyectiva(dim=2)/subespacio(dim=1,dim_ambiente=2) - eje de la involución (por defecto no se determina) \\
    # fijos: tupla(vector(3), vector(3)) - puntos fijos de la involución, ambos puntos de la cónica (por defecto no se determinan)
    #
    def __init__(self, conica, vertice = None, fijos = None, eje = None):
        assert vertice is not None or fijos is not None or eje is not None, "Se debe indicar alguna informacion para crear una involucion"
        self._es_conica = conica._tipo() == 0
        self._conica = conica
        if vertice is not None:
            self._vertice = vertice
            self._fijos = None
            self._eje = None
        else:
            if not self._es_conica:
                paso("Necesitaremos la conica sin parametrizar para lo siguiente:")
            c = conica if self._es_conica else conica.conica()
            if eje is not None:
                self._eje = eje
                self._fijos = None
            else:
                self._eje = None
                self._fijos = fijos
                paso("El vertice se calcula como el polo del eje:")
                e = subespacio(fijos[0], fijos[1])
                self._vertice = c.polo(e)
                paso("El vertice es: ", self._vertice)

    # Métodos accedentes

    #\m
    # Devuelve la cónica sobre la que se aplica esta involución. Puede ser de tipo conica o parametrizacion_conica.
    def conica(self):
        return self._conica

    #\m
    # Devuelve el vértice de esta involución.
    def vertice(self):
        return self._vertice

    # Otros métodos

    #\m
    # Devuelve el eje de esta involución.
    #
    # Implementación \\
    # Si se conocía en el momento de la creación simplemente lo devuelve. \\
    # Si no, calcula la polar del vértice.
    #
    def eje(self):
        if self._eje is not None:
            return self._eje
        paso("El eje es simplemente la polar del vertice")
        if self._es_conica:
            return self._conica.polar(self._vertice)
        paso("Recuperamos la conica sin parametrizar")
        return self._conica.conica().polar(self._vertice)

    #\m
    # Devuelve los dos puntos fijos de esta involución.
    #
    # Implementación \\
    # Si se conocían en el momento de la creación simplemente los devuelve. \\
    # Si no, calcula la intersección del eje con la cónica.
    #
    def puntos_fijos(self):
        if self._fijos is not None:
            return self._fijos
        paso("Los puntos fijos son la interseccion del eje con la conica")
        eje = self.eje()
        if self._es_conica:
            rep = eje.representantes()
            _no_pasos()
            r = recta_proyectiva(rep[0], rep[1])
            _no_pasos(False)
            theta = var('theta_var', latex_name = '\\theta')
            paso("Obtenemos una parametrizacion de la recta con: ", rep, "; ", r[theta])
            return self._conica.interseccion_recta(r)
        return self._conica.interseccion_recta(eje)

    #\m
    # Calcula la imagen mediante esta involución del punto dado. Si se tiene una parametrización de la cónica se puede
    # dar la coordenada no homogénea del punto.
    #
    # Uso: h(x) (donde h es una involución y x un punto de la cónica).
    #
    # Implementación \\
    # Calcula V(V, x) como subespacio o recta_proyectiva dependiendo de cómo sea haya dado la cónica y lo interseca con
    # esta. Se devuelve el otro punto de intersección.
    #
    # Parámetros \\
    # x: complejo/Infinity/vector(3) - punto de la cónica que se quiere calcular su imagen
    #
    def __call__(self, x):
        if es_parametro(x):
            assert not self._es_conica, "Solo se puede dar un parametro para calcular la imagen si la conica esta parametrizada"
            x = self._conica[x]
        else:
            assert x in self._conica, "El punto debe pertenecer a la conica"
        paso("Montamos la recta que pasa por el vertice y el punto dado para luego intersecar con la conica")
        r = recta_proyectiva(self._vertice, x) if self._es_conica else subespacio(self._vertice, x)
        (p, q) = self._conica.interseccion_recta(r)
        return p if matrix([p, x]).rank() == 2 else q

    def __repr__(self):
        return "<Involucion con vertice " + str(self._vertice) + " sobre la conica " + str(self._conica) + ">"

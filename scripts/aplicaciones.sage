#\s
# Archivo que contiene clases relativas a aplicaciones proyectivas con operaciones simples y funciones creadoras de estas.
#
# Autor: Pablo Sanz Sanz
#

# Funciones globales

#\f
# Crea un objeto del tipo aplicacion_proyectiva dadas las transformaciones que se quieren asegurar y su centro.
#
# Importante tener en cuenta que para los parámetros se piden listas de vectores, luego la llamada podría ser algo así: \\
# crear_aplicacion_proyectiva([A1, B1, C1], [A2, B2, C2], [Z]) con Ai, Bi, Ci, Z puntos de P^2 (3 coordenadas).
#
# Implementación \\
# Tomando como referencia los puntos iniciales junto con el centro se tiene que se transforman e_1, ..., e_k en cada uno de los
# puntos finales, y e_(k+1), ... e_(n+1) en 0 (como aplicación lineal). Por tanto, la matriz de la aplicación serán los puntos
# finales por columnas (multiplicados por los coeficientes que hacen que [f(e_1 + ... e_(k-1))] = f(e_k)) y cero para el resto de
# columnas. Después hay que deshacer el cambio de referencia multiplicando por la inversa de la matriz que forma la base elegida
# por la derecha.
#
# PROBABLEMENTE NO FUNCIONARÁ CON PUNTOS PARAMÉTRICOS.
#
# Parámetros \\
# iniciales: lista(vector(k<=n+1)) - puntos proyectivamente independientes que se quieren transformar \\
# finales: lista(vector(k<=n+1)) - puntos a los que se quiere llevar los puntos iniciales (en el mismo orden) \\
# centro: lista(vector(n+2-k)) - puntos proectivamente independientes a los iniciales que pertenecerán al centro
# de la aplicación proyectiva, esto es, al núcleo de la aplicación lineal asociada (por defecto vacía)
#
def crear_aplicacion_proyectiva(iniciales, finales, centro = []):
    assert len(iniciales) > 0, "Se debe indicar que puntos se quiere transformar"
    assert len(iniciales) == len(finales), "Para cada punto inicial se debe especificar su imagen"
    d = len(iniciales[0])
    assert len(iniciales) + len(centro) == d + 1, "Se deben especificar suficientes puntos entre inciales y centro para formar una referencia"
    assert all(map(lambda p: len(p) == d, iniciales + finales + centro)), "Todos los puntos deben pertenecer al mismo espacio"
    paso("Para crear la aplicacion vamos a usar los puntos inciales mas el centro como referencia (con el ultimo inicial como unidad)")
    paso("Calculamos la matriz del cambio de referencia a la canonica:")
    # Matriz del cambio de referencia a la canónica
    p = matriz_asociada(matrix(iniciales[0 : len(iniciales) - 1] + centro + [iniciales[-1]]).T)
    paso(p)
    coef = vector([var('alpha_var' + str(i + 1), latex_name = '\\alpha_' + str(i + 1)) for i in range(len(iniciales) - 1)])
    ec = (matrix(finales[0 : len(finales) - 1]).T * coef).simplify_full()
    paso("Forzamos [f(e_1 + ... + e_(k-1))] = f(e_k):")
    paso(matrix([ec]).T, "=", matrix([finales[-1]]).T)
    res = solve((ec - finales[-1]).simplify_full().list(), coef.list())
    # Finales * Alfa y rellenamos los demás con 0
    m = matrix([finales[i] * coef[i].substitute(res[0][i]) for i in range(len(finales) - 1)] \
                + [vector([0 for j in range(d)]) for i in range(d - len(finales) + 1)]).T
    paso("En la referencia formada por", iniciales[0 : len(iniciales) - 1] + centro + [iniciales[-1]], ", la matriz de la aplicacion es:")
    paso(m)
    paso("El cambio de referencia se hace:")
    paso(m, p, "^-1 = ", m * p^-1)
    return aplicacion_proyectiva(m * p^-1)

# Clases

#\c
# Clase que representa una aplicación proyectiva de un subespacio arbitario en otro.
#
class aplicacion_proyectiva:

    #\i
    # Construye una aplicacion poryectiva dadas su matriz asociada.
    #
    # En general, no está pensado para usarse directamente, sino para ser creada por otras funciones.
    #
    # Parámetros \\
    # matriz: matriz(n, m) - matriz asociada a la aplicacion
    #
    def __init__(self, matriz):
        self._matriz = matriz

    # Métodos accedentes

    #\m
    # Devuelve la matriz asociada a esta aplicación.
    def matriz_asociada(self):
        return self._matriz

    # Otros métodos

    #\m
    # Devuelve el centro de esta aplicación proyectiva, que es el núcleo de la lineal asociada.
    #
    # Implementación \\
    # Calcula el subespacio dual al que contiene las filas de la matriz asociada.
    #
    def centro(self):
        paso("El centro es el nucleo de la aplicacion lineal asociada")
        paso("Se obtiene directamente de las filas de la matriz asociada", self._matriz, "por filas, vistas como ecuaciones (dual)")
        _no_pasos()
        dual = subespacio(filter(lambda x: x != 0, self._matriz.rows()))
        _no_pasos(False)
        return dual.dual()

    #\m
    # Determina si esta aplicación es inyectiva.
    #
    # Implementación \\
    # Comprueba que el rango de la matriz asociada coincida con el núemro de columnas.
    #
    def es_inyectiva(self):
        return self._matriz.rank() == self._matriz.ncols()

    #\m
    # Determina si esta aplicación es sobreyectiva.
    #
    # Implementación \\
    # Comprueba que el rango de la matriz asociada coincida con el núemro de filas.
    #
    def es_sobreyectiva(self):
        return self._matriz.rank() == self._matriz.nrows()

    #\m
    # Determina si esta aplicación es una homografía (biyectiva).
    #
    # Implementación \\
    # Comprueba que sea inyectiva y la matriz cuadrada.
    #
    def es_homografia(self):
        return self._matriz.is_square() and self.es_inyectiva()

    #\m
    # Calcula la imagen mediante esta aplicación del punto dado.
    #
    # Uso: f(x) (donde f es una aplicacion_proyectiva, x un punto).
    #
    # Implementación \\
    # Multiplica el vector por la matriz asociada.
    #
    # Parámetros \\
    # x: vector(n) - punto del que se quiere calcular su imagen
    #
    def __call__(self, x):
        assert len(x) == self._matriz.ncols(), "El punto debe pertenecer al espacio inicial"
        _no_pasos()
        assert x not in self.centro(), "El punto no puede pertenecer al centro de la aplicacion"
        _no_pasos(False)
        return self._matriz * x

    #\m
    # Calcula los autovalores de esta aplicación, asumiendo que va de un espacio en sí mismo.
    #
    # Implementación \\
    # Resuelve det(f - lambda*id) = 0 para lambda (donde f es esta aplicación).
    #
    def autovalores(self):
        assert self._matriz.ncols() == self._matriz.nrows(), \
                "Para calcular autovalores la aplicacion debe ser de un espacio en si mismo"
        lambda0 = var('lambda_var', latex_name = '\\lambda')
        paso("det", self._matriz - lambda0, "=0")
        return map(lambda sol: sol.rhs(), solve((self._matriz - lambda0).det(), lambda0))

    #\m
    # Devuelve tuplas con los autovectores de la matriz asociada a esta aplicación. \\
    # Para cada tupla, su primer elemento es el autovalor, el segundo la base del subespacio y el tercero la multiplicidad. \\
    #
    # Implementación \\
    # Utilzia eigenvectors_right() de Sage. Para obtener los puntos fijos solamente usar puntos_fijos().
    #
    def autovectores(self):
        return self._matriz.eigenvectors_right()

    #\m
    # Devuelve los puntos fijos de esta aplicacion proyectiva.
    #
    # Implementación \\
    # Para cada autovalor lambda resuelve (f - lambda id)X = 0 para X (donde f es esta aplicación).
    #
    def puntos_fijos(self):
        _no_pasos()
        autovalores = self.autovalores()
        _no_pasos(False)
        lambda0 = var('lambda_var', latex_name = '\\lambda')
        paso("Resolvemos", self._matriz - lambda0, "*X = 0 para los autovalores:", self.autovalores())
        vars = [var('x_var' + str(i), latex_name = 'x_' + str(i)) for i in range(self._matriz.ncols())]
        x = vector(vars)
        paso("Si queda algun parametro, se podra sacar como factor comun y sustituirse por cualquier valor")
        # Resuelve el sistema y se queda sólo con los resultados como vectores
        sol = map(lambda sol: vector(map(lambda x: x.rhs(), sol)).simplify_full(), \
                [solve(((self._matriz - autovalor) * x).list(), vars)[0] for autovalor in autovalores])
        # Divide entre el MCD por si acaso
        return map(lambda p: p / gcd(p.list()), sol)
        #return map(lambda sol: vector(map(lambda x: x.rhs(), sol[0])).simplify_full(), \
        #        [solve(((self._matriz - autovalor) * x).list(), vars) for autovalor in autovalores])

    #\m
    # Devuelve una nueva aplicación con los cambios de referencia especificados. Si R1 y R2 son las referencias de cada
    # uno de los espacios y R1', R2' las nuevas, se pasarán como parámetro las matrices del cambio de R1 a R1' y de R2 a R2'.
    #
    # Implementación \\
    # Si las matrices son M, P y Q las de esta aplicación, el cambio de coordenadas incial y final, respectivamente,
    # calcula Q*M*P^-1.
    #
    # Parámetros \\
    # cambio_inicial: matriz(m, m) - matriz del cambio de referencia del espacio de salida (por defecto I) \\
    # cambio_final: matriz(n, n) - matriz del cambio de referencia del espacio de llegada
    # (por defecto I o la misma que la anterior) \\
    # misma: booleano - determina si se debe usar el mismo cambio de base en ambos espacios en caso de que no se especifique
    # la matriz del espacio de llegada y sea posible
    #
    def cambiar_referencias(self, cambio_inicial = None, cambio_final = None, misma = True):
        if cambio_inicial is None:
            cambio_inicial = 1
        if cambio_final is None:
            if misma and self._matriz.is_square():
                cambio_final = cambio_inicial
            else:
                cambio_final = 1
        return aplicacion_proyectiva(cambio_final * self._matriz * cambio_inicial^-1)

    #\m
    # Operador *. Devuelve la composición de las aplicaciones ((self o otra)(x) = self(otra(x))).
    #
    # Uso: f * g (f y g son aplicaciones compatibles para la composición).
    #
    # Implementación \\
    # Devuelve una nueva aplicación cuya matriz sea el producto de las de ambas.
    #
    # Parámetros \\
    # otra: aplicacion_proyectiva - aplicación con la que componer
    #
    def __mul__(self, otra):
        assert self._matriz.ncols() == otra._matriz.nrows(), \
                "El espacio de llegada de la primera aplicacion debe ser el de salida de la segunda"
        return aplicacion_proyectiva(self._matriz * otra._matriz)

    #\m
    # Operador ^ (ó **). Devuelve el resultado de componer una aplicación consigo misma n veces (^-1 devuelve la inversa).
    #
    # Uso f^n (ó f**n) (f es una homografía de la recta y n un entero).
    #
    # Implementación \\
    # Devuelve una nueva aplicación cuya matriz es la de esta elevada a n.
    #
    # Parámetros \\
    # n: entero - exponente al que elevar
    #
    def __pow__(self, n):
        assert self._matriz.is_square(), "La aplicacion debe ir de un espacio en si mismo"
        return aplicacion_proyectiva(self._matriz^n)

    def __repr__(self):
        return "<Aplicacion proyectiva con matriz asociada\n" + str(self._matriz) + ">"

#\c
# Clase que representa una proyección dados un subespacio centro y un subespacio imagen.
#
class proyeccion:

    #\i
    # Construye una proyección dados su centro (Z) e imagen (Y) de un mismo espacio Z.
    # Se debe asegurar que dim Z + dim Y = dim X - 1.
    #
    # Parámetros \\
    # centro: subespacio - centro de la proyección \\
    # imagen: subespacio - espacio de llegada de la aplicación
    #
    def __init__(self, centro, imagen):
        dim = centro.dimension_ambiente()
        assert dim == imagen.dimension_ambiente(), "Los subespacios deben pertenecer al mismo espacio ambiente"
        assert centro.dim() + imagen.dim() == dim - 1,  "Se debe cumplir dim Z + dim Y = dim X - 1"
        self._centro = centro
        self._imagen = imagen

    # Métodos accedentes

    #\m
    # Devuelve el centro de esta proyección.
    def centro(self):
        return self._centro

    #\m
    # Devuelve el espacio de llegada de esta proyección.
    def imagen(self):
        return self._imagen

    # Otros métodos

    #\m
    # Calcula la imagen mediante esta proyección del punto dado.
    #
    # Uso: pi(x) (donde pi es una proyeccion y x un punto).
    #
    # Implementación \\
    # Calcula el subespacio V(Z, x) y lo interseca con el subespacio de llegada, asumiendo x no pertenece a Z.
    #
    # Parámetros \\
    # x: vector(n) - punto del que se quiere calcular su imagen
    #
    def __call__(self, x):
        _no_pasos()
        assert x not in self._centro, "Los puntos del centro no tienen imagen en una proyeccion"
        punto = subespacio(x)
        _no_pasos(False)
        paso("Calculamos V(Z, x)")
        V = self._centro.suma(punto)
        paso("Ahora intersecamos con el espacio de llegada, lo que nos deberia dar un solo punto")
        # La intersección debería ser un punto
        return V.interseccion(self._imagen).representantes()[0]

    #\m
    # Devuelve esta proyección como un objeto del tipo aplicacion_proyectiva.
    #
    # Implementación \\
    # Utiliza la función creadora de aplicaciones proyectivas crear_aplicacion_proyectiva, dando como puntos transformados
    # los representantes del subespacio imagen (que se transforman en sí mismos) y como centro los representantes del centro
    # de esta proyección. Como punto unidad escoge la suma de los representantes del espacio imagen y el centro y calculamos
    # su imagen.
    #
    def aplicacion(self):
        _no_pasos()
        repr_imagen = self._imagen.representantes()
        repr_centro = self._centro.representantes()
        unidad = reduce(lambda p1, p2: p1 + p2, repr_imagen) + reduce(lambda p1, p2: p1 + p2, repr_centro)
        pi_unidad = self(unidad)
        _no_pasos(False)
        paso("Los puntos: ", repr_imagen, " se transforman en si mismos y ", repr_centro, " son el centro de la aplicacion")
        paso("El punto unidad:", unidad, " se transforma en ", pi_unidad)
        return crear_aplicacion_proyectiva(repr_imagen + [unidad], repr_imagen + [pi_unidad], repr_centro)

    def __repr__(self):
        return "<Proyeccion de centro " + str(self._centro) + " e imagen " + str(self._imagen) + ">"

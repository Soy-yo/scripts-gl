# \s
# Este archivo proporciona diferentes funciones y clases para trabajar con espacios proyectivos de forma sencilla.
#
# Las implementaciones no son eficientes y seguro que se pueden mejorar, pero, en general, los vectores que se van a usar no van
# a tener una dimensión muy grande, por lo que no debería notarse demasiado la lentitud.
#
# Autor: Pablo Sanz Sanz
#

# Funciones globales

#\f
# Devuelve la matriz del cambio de referencia entre la inicial y la final.
#
# Implementación \\
# Calcula las matrices asociadas a cada referencia (ini y fin) y devuelve fin^-1 * ini
# (pasamos la incial a la canónica y de esta a la final). ¡Recordar multiplicar el vector
# por la derecha!
#
# Parámetros \\
# inicial: matriz(n, n+1) - matriz representante del sistema de referencia original \\
# final: matriz(n, n+1) - matriz representante del sistema de referencia de destino
#
def cambiar_referencia(inicial, final):
    assert inicial.dimensions() == final.dimensions(), \
           "Las referencias deben pertenecer al mismo espacio proyectivo"
    ini = matriz_asociada(inicial)
    fin = matriz_asociada(final)
    paso("Vamos a pasar primero de la inicial a las coordenadas canonicas y despues a la final.")
    paso("Calculamos B^-1 * A (B y A asociadas a final e inicial, respectivamente):")
    paso(fin^-1, ini, "=", fin^-1 * ini)
    return fin^-1 * ini

#\f
# Devuelve una matriz representando una referencia canónica con punto unidad la suma de todos ellos.
#
# Pensado para uasrse al cambiar de referencia.
#
# Parámetros \\
# n: entero - dimensión del espacio (la matriz será de (n+1)x(n+2))
#
def referencia_canonica(n):
    return identity_matrix(n + 1).insert_row(n + 1, vector([1 for i in range(n + 1)])).T

#\f
# Devuelve la matriz asociada a una referencia dada por una matriz n x n+1.
#
# Implementación \\
# Monta el sistema a_0*x_0 + ... + a_n*x_n = x_(n+1), donde x_i representa el punto
# i-ésimo de la referencia, y devuelve la misma matriz eliminando la última columna y
# multiplicada cada columna por a_i.
#
# Parámetros \\
# ref: matriz(n, n+1) - matriz representante del sistema de referencia
#
def matriz_asociada(ref):
    assert es_referencia(ref), "La matriz debe ser una referencia (nx(n+1) y rango maximo)"
    m = ref.matrix_from_columns(range(ref.ncols() - 1))
    coef = vector([var('alfa' + str(i + 1), latex_name = r'alfa_' + str(i + 1)) for i in range(0, m.nrows())])
    ec = m * coef
    paso("Forzamos [x_0 + ... + x_n] = x_{n+1}:")
    paso(matrix([ec]).T, "=", matrix([ref.columns()[-1]]).T)
    res = solve((ec - ref.columns()[-1]).list(), coef.list())
    return matrix([m.columns()[i] * coef[i].substitute(res[0][i]) for i in range(0, m.ncols())]).T

#\f
# Determina si la matriz dada puede representar una referencia proyectiva.
#
# Implementación \\
# Comprueba que sea una matriz cuadrada más la columna del punto unidad y que el rango sea máximo.
#
def es_referencia(matriz):
    return matriz.ncols() == matriz.nrows() + 1 and \
            matriz.matrix_from_columns(range(matriz.ncols() - 1)).det() != 0

# Clases

#\c
# Clase que representa un subespacio arbitrario en un espacio de dimensión arbitraria.
# Las dimensiones son determinadas por la dimensión de los puntos dados y por la independencia de los mismos.
#
# Si se quiere crear un subespacio a partir de unas ecuaciones hay que extraer los coeficientes de estas como vectores, crear
# un subespacio con ellos y acceder a su dual. Es una operación bastante más costosa, pero ahorra tiempo de programación.
# Ejemplo. Recta 2x-y+1=0, proyectivamente 2x-y+z=0, se representa en el dual como (2 : -1 : 1). Así se instanciaría como
# subespacio(vector([2, -1, 1])).dual().
#
class subespacio:

    #\i
    # Inicializa los atributos.
    #
    # Implementación \\
    # Se calculan la dimensión ambiente y la del subespacio. Esta última mediante el rango (menos 1) de la matirz formada por
    # los puntos dados, una vez reducida la matriz a su forma de Hermite. De esta matriz, multiplicada por unos parámetros
    # lambda(i) se obtienen las ecuaciones paramétricas. Para las implícitas se completa esta matriz con las coordenadas x(i)
    # y se igualan a 0 los determinantes de las submatrices formadas por las primeras dim+1 filas más una fila extra hasta
    # lograr una matriz cuadrada (tomadas por orden).
    #
    # Parámetros \\
    # *puntos: vector(n+1) - puntos que contiene el subespacio generado
    #
    def __init__(self, *puntos):
        if len(puntos) == 0 or len(puntos[0]) == 0:
            self.__subespacio_vacio()
        else:
            # Si recibe una lista en vez de varargs
            if type(puntos[0]) is list:
                puntos = puntos[0]
            self._dimension_ambiente = self.__dimension_ambiente(puntos)
            assert 0 not in puntos, "0 no es un punto proyectivo"
            assert self._dimension_ambiente >= 0, "Todos los puntos deben pertenecer al mismo espacio"
            m = matrix([p.list() for p in puntos]).echelon_form().T
            paso("Subespacio que contiene los puntos ", puntos)
            paso("Primero reducimos a su forma normal: ", m)
            self._dim = m.rank() - 1
            self._matriz = m.matrix_from_columns(range(0, self._dim + 1))
            self._params = self.__params()
            self._vars = self.__vars()
            self._parametricas = self.__parametricas()
            self._implicitas = self.__implicitas()

    # Métodos accedentes

    #\m
    # Dimensión del espacio ambiente de este subespacio.
    def dimension_ambiente(self):
        return self._dimension_ambiente

    #\m
    # Dimensión de este subespacio.
    def dim(self):
        return self._dim

    #\m
    # Ecuaciones paramétricas de este subespacio.
    def parametricas(self):
        return self._parametricas

    #\m
    # Devuelve un vector conteniendo los parámetros de la expresión paramétrica de este subespacio.
    def parametros(self):
        return self._params

    #\m
    # Ecuaciones implícitas de este subespacio.
    def implicitas(self):
        return self._implicitas

    #\m
    # Devuelve las variables de las ecuaciones implícitas de este subespacio.
    def variables(self):
        return self._vars;

    #\m
    # Devuelve dim + 1 vectores independientes que pertenecen al subespacio.
    def representantes(self):
        return self._matriz.columns()

    #\m
    # Si este subespacio es un único punto devuelve este punto como objeto vector (equivalente a self.representantes()[0])
    def punto(self):
        assert self.es_punto(), "El subespacio debe ser un unico punto"
        return self.representantes()[0]

    #\m
    # Determina si este subespacio es vacío.
    def es_vacio(self):
        return self._dimension_ambiente == -1

    #\m
    # Determina si este subespacio es un único punto.
    def es_punto(self):
        return self._dim == 0

    #\m
    # Determina si el subespacio coincide con su espacio ambiente.
    def es_total(self):
        return self._dim == self._dimension_ambiente

    # Otros métodos

    #\m
    # Devuelve el subespacio dual a este subespacio.
    #
    # Implementación \\
    # Obtiene los coeficientes de las ecuaciones implícitas y los usa como vectores para el nuevo subespacio.
    #
    def dual(self):
        return subespacio([self.__vector_dual(ec) for ec in self._implicitas])

    #\m
    # Devuelve V(self, otro).
    #
    # Implementación \\
    # Crea un subespacio uniendo los representantes de los dados.
    #
    # Parámetros \\
    # otro: subespacio - subespacio con el que sumar este
    #
    def suma(self, otro):
        assert self.__mismo_ambiente(otro), "Los subespacios a sumar deben pertenecer al mismo espacio"
        paso("Para unir creamos un subespacio con los representantes de cada uno: ", [self.representantes(), otro.representantes()])
        return subespacio(self.representantes() + otro.representantes())

    #\m
    # Devuelve self (intersección) otro
    #
    # Implementación \\
    # Crea un subespacio uniendo los representantes de los duales dados (devuelve el dual).
    #
    # Parámetros \\
    # otro: subespacio - subespacio con el que intersecar este
    #
    def interseccion(self, otro):
        assert self.__mismo_ambiente(otro), "Los subespacios a intersecar deben pertenecer al mismo espacio"
        if self.es_vacio() or otro.es_vacio():
            return subespacio()
        paso("Para intersecar unimos las ecuaciones de cada subespacio: ", [self.implicitas(), otro.implicitas()])
        _no_pasos()
        dual1 = self.dual()
        dual2 = otro.dual()
        _no_pasos(False)
        s = subespacio(dual1.representantes() + dual2.representantes())
        paso("Ahora obtenemos el dual del subespacio calculado (era el dual del que buscabamos):")
        res = s.dual()
        return res

    #\m
    # Devuelve las ecuaciones inhomogéneas del subespacio con xn == 1.
    #
    # Implementación \\
    # Sustituye xn == 1 en todas las ecuaciones.
    #
    def ecuaciones_inhomogeneas(self):
        return [ec.substitute(self._vars[-1] == 1) for ec in self._implicitas]

    #\m
    # Operador in. Determina si un punto está contenido en este subespacio o no.
    #
    # Uso: P in U (P es un punto y U un subespacio).
    #
    # Implementación \\
    # Sustituye las coordenadas del punto en cada una de las ecuaciones implícitas y devuelve si todas se complen.
    #
    # Parámetros \\
    # punto: vector(n) - punto que comprobar si pertenece al subespacio
    #
    def __contains__(self, punto):
        if self.es_vacio():
            return False
        assert len(punto) == self._dimension_ambiente + 1, "El punto debe pertenecer al espacio ambiente"
        subs = [self._vars[i] == punto[i] for i in range(len(punto))]
        ecs = [ec.substitute(subs) for ec in self._implicitas]
        return solve(ecs, self._vars.list())

    #\m
    # Operador ==. Determina si dos subespacios son iguales.
    #
    # Uso: U == V (U y V son subespacios).
    #
    # Implementación \\
    # Comprueba que las dimensiones coincidan y que al sumarlos quede un subespacio con la misma dimensión (se usa la suma
    # porque la intersección es más lenta).
    #
    # Parámetros \\
    # otro: subespacio - subespacio a comprobar la igualdad
    #
    def __eq__(self, otro):
        assert self.__mismo_ambiente(otro), "Los subespacios a comprobar la igualdad deben pertenecer al mismo espacio"
        # La dimensión debe ser la misma y al sumarlos el subespacio no debería cambiar
        _no_pasos()
        res = self._dim == otro._dim and self.suma(otro)._dim == self._dim
        _no_pasos(False)
        return res

    #\m
    # Operador !=. Determina si dos subespacios son diferentes.
    #
    # Uso: U != V (U y V son subespacios).
    #
    # Implementación \\
    # Negar el operador ==.
    #
    # Parámetros \\
    # otro: subespacio - subespacio a comprobar la desigualdad
    #
    def __ne__(self, otro):
        return not (self == otro)

    # Métodos auxiliares

    def __subespacio_vacio(self):
        self._dimension_ambiente = -1
        self._dim = -1
        self._matriz = matrix()
        self._params = []
        self._vars = []
        self._parametricas = []
        self._implicitas = []

    def __dimension_ambiente(self, puntos):
        if len(puntos) == 0:
            return -1
        dim = len(puntos[0].list())
        # Todos los puntos deben tener el mismo tamaño
        for p in puntos:
            if len(p.list()) != dim:
                return -1
        return dim - 1

    def __mismo_ambiente(self, otro):
        return self.es_vacio() or otro.es_vacio() or self._dimension_ambiente == otro._dimension_ambiente

    def __params(self):
        return vector([var('lambda' + str(i + 1), latex_name = r'lambda_' + str(i + 1)) for i in range(self._dim + 1)])

    def __vars(self):
        return vector([var('x' + str(i + 1), latex_name = 'x_' + str(i + 1)) for i in range(self._dimension_ambiente + 1)])

    def __parametricas(self):
        arr = (self._matriz * self._params).list()
        return [self._vars[i] == arr[i] for i in range(len(arr))]

    def __implicitas(self):
        paso("Resolvemos las siguientes ecuaciones:")
        arr = self.__mover_ceros(self.__unir_vars(self.__array_2d(self._matriz)))
        ec = []
        for i in range(self._dim + 1, self._dimension_ambiente + 1):
            m = matrix(arr[0 : self._dim + 1] + [arr[i]])
            ec.append(m.det().factor().simplify_full() == 0)
            paso("det", m, "= 0 => ", ec[-1])
        return ec

    def __array_2d(self, matriz):
        return [f.list() for f in matriz.rows()]

    def __unir_vars(self, array):
        for i in range(len(array)):
            array[i].append(self._vars[i])
        return array

    def __mover_ceros(self, array):
        i = 0
        j = self._dim + 1
        while i < self._dim + 1:
            if vector(array[i][0 : -1]) == 0:
                temp = array[j]
                array[j] = array[i]
                array[i] = temp
                j = j + 1
            else:
                i = i + 1
        return array

    def __vector_dual(self, ecuacion):
        # Nueva implementación: es más clara y creo que da menos lugar a errores
        v = [ecuacion.lhs().substitute(map(lambda x: x == 0, filter(lambda xj: bool(xj != xi), self._vars)), xi == 1) for xi in self._vars]
        return vector(v)
        '''
        temp = []
        # Variable inventada para ecuaciones del tipo a*x_i == 0 que, tras sustituir devuelve un operands() vacío
        var('x_artificial')
        for xi in self._vars:
            if xi in ecuacion.free_variables():
                # Aprovecha que, tras sustituir, el término independiente va al final
                temp.append((ecuacion.substitute(xi == 1).lhs() + x_artificial).operands()[-1])
            else:
                temp.append(0)
        return vector(temp)
        '''

    def __repr__(self):
        if self.es_vacio():
            return "<Subespacio vacio>"
        if self.dim() == 0:
            return "<Punto " + str(self.representantes()[0]) + ">"
        if self.es_total():
            return "<Espacio de dimension " + str(self._dim) + ">"
        return "<Subespacio de dimension " + str(self._dim) + " en un espacio de dimension " + str(self._dimension_ambiente) \
                + " con ecuacion(es) implicita(s) " + str(self._implicitas) + ">"

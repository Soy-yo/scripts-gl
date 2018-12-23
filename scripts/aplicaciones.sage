#\s
# Archivo que contiene clases relativas a aplicaciones proyectivas con operaciones simples y funciones creadoras de estas.
#
# Autor: Pablo Sanz Sanz
#

# Clases

#\c
# Clase que representa una aplciación proyectiva de un subespacio arbitario en otro.
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
        ocultar_procedimiento()
        dual = subespacio(self._matriz.rows())
        reanudar_procedimiento()
        return dual.dual()
    
    #\m
    # Determina si esta apliación es inyectiva.
    #
    # Implementación \\
    # Comprueba que el rango de la matriz asociada coincida con el núemro de columnas.
    #
    def es_inyectiva(self):
        return self._matriz.rank() == self._matriz.ncols()
    
    #\m
    # Determina si esta apliación es sobreyectiva.
    #
    # Implementación \\
    # Comprueba que el rango de la matriz asociada coincida con el núemro de filas.
    #
    def es_sobreyectiva(self):
        return self._matriz.rank() == self._matriz.nrows()
    
    #\m
    # Determina si esta apliación es una homografía (biyectiva).
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
        ocultar_procedimiento()
        assert x not in self.centro(), "El punto no puede pertenecer al centro de la aplicacion"
        mostrar_procedimiento()
        return self._matriz * x
        
    #\m
    # Calcula los autovalores de esta aplicación, asumiendo que va de un espacio en sí mismo.
    #
    # Implementación \\
    # Resuelve det(f - lambda*id) = 0 para lambda (donde f es esta apliación).
    #
    def autovalores(self):
        assert self._matriz.ncols() == self._matriz.nrows(), \
                "Para calcular autovalores la aplicacion debe ser de un espacio en si mismo"
        var('lambda0', latex_name = r'lambda')
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
    # Para cada autovalor lambda resuelve (f - lambda id)X = 0 para X (donde f es esta apliación).
    #
    def puntos_fijos(self):
        ocultar_procedimiento()
        autovalores = self.autovalores()
        reanudar_procedimiento()
        var('lambda0', latex_name = r'lambda')
        paso("Resolvemos", self._matriz - lambda0, "*X = 0 para los autovalores:", self.autovalores())
        vars = [var('x' + str(i)) for i in range(self._matriz.ncols())]
        x = vector(vars)
        paso("Si queda algun parametro, se podra sacar como factor comun y sustituirse por cualquier valor")
        # Resuelve el sistema y se queda sólo con los resultados como vectores
        return map(lambda sol: vector(map(lambda x: x.rhs(), sol[0])).simplify_full(), \
                [solve(((self._matriz - autovalor) * x).list(), vars) for autovalor in autovalores])
    
    #\m
    # Devuelve una nueva apliación con los cambios de referencia especificados.
    #
    # Implementación \\
    # REVISAR. Si las matrices son M, P y Q las de esta aplciación, el cambio de coordenadas incial y final, respectivamente,
    # calcula Q*M*P^-1. Posiblemente esté al revés, esto es lo que hay que revisar.
    #
    # Parámetros \\
    # cambio_inicial: matriz(m, m) - matriz del cambio de referencia del espacio de salida (por defecto I) \\
    # cambio_final: matriz(n, n) - matriz del cambio de referencia del espacio de llegada
    # (por defecto I o la misma que la anterior) \\
    # misma: booleano - determina si se debe usar el mismo cambio de base en ambos espacios en caso de que no se especifique
    # la matriz del espacio de llegada y sea posible
    #
    def cambiar_referencias(self, cambio_inicial = None, cambio_final = None, misma = True):
        if cambio_inicial == None:
            cambio_inicial = 1
        if cambio_final == None:
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
    # otra: aplicacion_proyectiva - aplciación con la que componer
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
    # Devuelve una nueva aplciación cuya matriz es la de esta elevada a n.
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
    # centro: subespacio - centro de la proyección
    # imagen: subespacio - espacio de llegada de la apliación
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
        ocultar_procedimiento()
        assert x not in self._centro, "Los puntos del centro no tienen imagen en una proyeccion"
        punto = subespacio(x)
        reanudar_procedimiento()
        paso("Calculamos V(Z, x)")
        V = self._centro.suma(punto)
        paso("Ahora intersecamos con el espacio de llegada, lo que nos deberia dar un solo punto")
        # La intersección debería ser un punto
        return V.interseccion(self._imagen).representantes()[0]
        
    def __repr__(self):
        return "<Proyeccion de centro " + str(self._centro) + " e imagen " + str(self._imagen) + ">"
        
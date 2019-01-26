# scripts-gl
Scripts que pueden resultar útiles para resolver ejercicios de Geometría Lineal usando Sage.

Están pensados únicamente para su uso en resolución de ejercicios y exámenes de la asignatura GL de la Facultad de Matemáticas de la Universidad Complutense de Madrid. Por ello, no se pretende que los algoritmos sean óptimos ni 100% fiables, aunque sí es recomendable que tengan una fiabilidad alta.

Cualquier error que se detecte debería ser comunicado lo antes posible.

# Importante
Todos estos scripts han sido programados y probados por mí mismo en apenas unos meses. He intentado probar que todo funcione correctamente y he ido arreglando todos los errores que he ido encontrando. Parece que todo funciona correctamente. Pero eso no significa que todo funcione correctamente. Así que <b>no me hago responsable de posibles errores que puedan suceder a la hora de usar los scripts en un examen</b>. Se recomienda comprobar que las soluciones tienen sentido o incluso dibujarlas con GeoGebra.

# Fecha de salida
Salvo arreglos de última hora, que serán avisados, se asegura que a partir del <b>lunes 28 de enero a las 10:00</b> ya no habrá más modificaciones de código.

# Cómo usarlos
Las carpetas importantes son scripts y doc. En scripts se encuentran todos los archivos con scripts disponibles. En doc se encuentra la documentación básica para ser capaz de usar dichos scripts, sin necesidad de mirar ningún tipo de código.

<ol>
  <li>
    Descargamos el proyecto (o tan sólo las carpetas <b>scripts</b> y <b>doc</b> y recomendable también este <b>README</b>).
  </li>
  <li>
    En el examen los llevamos en un pendrive y los copiamos al ordenador. Abrimos el cuaderno. Si reconoce la carpeta directamente nos podemos ahorrar esto, si no, depende qué versión usemos:
  </li>
  <li>
    <ul>
      <li>
        En Jupyter (icono naranja) cargamos los archivos .sage con "Upload", creamos una nueva carpeta llamada "scripts" y los movemos ahí. Después creamos una nueva hoja en la que trabajar (en el directorio justo anterior a scripts (es decir, que podamos ver la carpeta y el archivo a la vez)) y ejecutamos en la primera línea <code>load('scripts/scripts_gl.sage')</code>. Con esto se cargarán todos los archivos y ya podremos trabajar. <i>Procedimiento recomendado</i>
      </li>
      <li>
        <strike>En la versión antigua (icono azul (o blanco?)) creamos directamente una nueva hoja y copiamos, pegamos y ejecutamos el contenido de cada uno de los archivos en las primeras filas. Después, ya se podrán usar sin más.</strike>
        Ya no están pensados para la versión antigua porque algunos archivos se enlazan mutuamente. Igualmente, se pueden unir todos los archivos en un mismo archivo para luego copiar y pegar en nuestro cuaderno. Esto debería funcionar, pero probablemente ralentice el ordenador.
      </li>
    </ul>
  </li>
  <li>
    Para ver los resultados con tipografía de Latex se recomienda ejecutar la línea <code>pretty_print_default(True)</code>.
  </li>
  <li>
    Abrimos <b>doc/index.html</b> en el navegador para tener toda la documentación necesaria.
  </li>
</ol>

Se ha añadido una nueva funcionalidad que permite mostrar por pantalla los pasos que va siguiendo una función, de forma que no sólo obtengamos el resultado sin saber cómo se ha calculado eso. Para más información ver la sección <b>Procedimientos</b> más abajo.

# Cómo participar
Si quieres subir algún script a este banco de scripts, puedes hacerlo sencillamente con un Pull-request a una rama nueva, o enviármelo, con unas ciertas condiciones.
<ol>
<li>
  Obviamente, el script debe funcionar y, si hay algún caso concreto que falla, se debe indicar en su documentación.
</li>
<li>
  Debe estar escrito como una función (o una clase con métodos) de Python. Esto es, debe recibir unos parámetros y recibir un resultado directamente.
  
<code>

    def foo(bar):
        ~ interior de la función ~
    
</code>
  
  No vale tener un cuaderno en el que sustituyes unos valores para unas variables y se va ejecutando paso a paso. De todas formas no está de más añadir, además, algo semejante si se quiere tener resultados intermedios.
</li>
<li>
  Como Python es como es, es obligatorio indentar, y toda la indentación debe ser uniforme. Así que se indentará con <b>4 espacios</b>.
</li>
<li>
  Cada función o clase debe tener justo sobre ella un comentario <code># bla bla</code> con un formato que se explciará más abajo explicando lo que hace y cómo lo hace. Esto es importante, pues en un examen se nos va a pedir algún tipo de procedimiento y debemos saber qué es lo que está haciendo la función por dentro. También debe indicar qué significa cada parámetro que recibe.
  
<code>

    # \f
    # Implentación \\
    # La función funciona así.
    #
    # Parámetros \\
    # a: vector(n) - vectorcito de entrada
    # b: recta - recta horrible con la que haremos cosas
    #
    def foo(a, b):
    
</code>

</li>
<li>
  Por simplicidad para todo el mundo, tanto el código como los comentarios deben estar todos en español (a pesar de que a mí, personalmente, me haga daño a la vista ver código en español).
</li>
</ol>

Si alguien quiere participar, pero no se ve capaz de implementar algo así, siempre me puede pedir ayuda para hacerlo.

# Procedimientos
Un tema importante de cara a un examen es el procedimiento seguido al resolver un ejercicio, pues de nada sirve poner directamente el resultado. Ante esto, he añadido la opción de mostrar los pasos que van dando las funciones, por ejemplo "Calculamos el determinante de M" (y añadir después el cálculo explícito). Esto ya depende de quién implemente cada función, si los añade o no. Tampoco es necesario especificar absolutamente todo lo que se hace porque luego no se entiende nada de la salida.

El caso, para poder ver los procedimientos hay que activarlo a mano: simplemente ejecutamos la siguiente instrucción: <code>mostrar_procedimiento(True)</code> y a partir de entonces las funciones comenzarán a mostrar la información que haya añadido el programador. Con esto activado hay que tener cuidado de no ejecutar muchas operaciones a la vez porque de nuevo se nos puede llenar la salida de operaciones intermedias y no entender nada. Por ejemplo, <code>subespacio(vector([1, 2, 1, 0]), vector([3, 4, 1, -1])).interseccion(subespacio(vector([1, 5, -2, 3])))</code> llena más de una pantalla de operaciones intermedias, pero podemos partirlo en varias operaciones y así será más claro.

Notemos que si no podemos cargar los archivos .sage debemos copiar y pegar el código del archivo procedimientos.sage en nuestro cuaderno.

Si somos programadores y consideramos que no es del todo obvio cómo se implementa nuestra función podemos utilizar la función <code>paso(*objetos)</code>. Su uso es muy sencillo: escribimos esta línea en el punto que queramos dar información, y como objetos pasamos todo lo que queramos mostrar, separados por comas: desde cadenas de texto o números hasta matrices o ecuaciones. Cualquier cosa. El resultado será una línea que se muestra por pantalla (solamente si el usuario activó previamente los procedimientos) con tipografía de Latex. Por ejemplo, <code>paso(A, B, "=", A * B)</code> muestra explícitamente la matriz A, la matriz B, el signo = y el resultado de su producto.

Importante no poner tildes ni caracteres raros porque es probable que el cuaderno no sea capaz de leerlo (comprobado con Jupyter 8.3).

# Ejemplos
Recientemente se ha añadido una nueva carpeta con cuadernos de Jupyter con algunos ejemplos de uso de los scripts. Para probarlos sólo hay que descargarlos y abrirlos con un cuaderno Jupyter. Se intentarán completar algunos ejercicios o exámenes utilizando únicamente los scripts.

# Formato de la documentación
En este proyecto se incluye un pequeño programa que se encarga de leer los comentarios de los archivos .sage para transformarlos en archivos .html (páginas web) que se puedan abrir con un navegador y sean más legibles que los propios comentarios en el código. Pero para que esto fucnione correctamente se necesite que los comentarios sigan unos estándares (sólo se refiere a los de documentación de las funciones; comentarios internos aclarativos pueden llevar el formato que sea).
<ul>
<li>
  Si se está documentando un archivo en general, es decir, lo que hay dentro de ese módulo, la primera línea del comentario debe ser \s. También debe haber una línea en blanco después de todo el comentario.
  
<code>
  
    # \s
    # Este archivo es muy bonito.
    # Hace cosas.
    #
    
    ~ código ~
    
</code>
  
</li>
<li>
  Para una función o método de una clase se utiliza \f y \m, respectivamente, aunque realmente no suponen ninguna diferencia real a la hora de generar la documentación.
  
<code>

    # \f
    # Esta función no hace nada.
    #
    def foo():
    
</code>
</li>
<li>
  Para una clase se utilizará \c
  
<code>
  
    # \c
    # Clase sin clase.
    #
    class clase:
    
</code>
  
</li>
<li>
  Para el constructor de una clase \i
  
<code>
  
    # \i
    # Constructor.
    #
    def __init__(self):
    
</code>
  
</li>
</ul>
Si en el archivo .html que se genere queremos que haya algún salto de línea en un lugar concreto tenemos dos posibilidades: dejar una línea entera en blanco (pero con comentario) o escribir \\. Más o menos como en Latex.

<code>

    # Esto se verá en la primera línea.
    # Esto también.
    # 
    # Pero esto ya en la segunda. \\
    # Y esto en la tercera.
  
</code>

import java.io.*
import java.lang.StringBuilder

const val TITLE = "Scripts de GL - Doc"
const val INDEX = "index.html"

const val NO_DOC_TYPE = 0
const val FILE_DOC = "\\s"
const val FILE_DOC_TYPE = 1
const val FUNCTION_DOC = "\\f"
const val FUNCTION_DOC_TYPE = 2
const val CLASS_DOC = "\\c"
const val CLASS_DOC_TYPE = 3
const val INIT_DOC = "\\i"
const val INIT_DOC_TYPE = 4
const val METHOD_DOC = "\\m"
const val METHOD_DOC_TYPE = 5

fun main() {
    val bw = BufferedWriter(FileWriter("doc/$INDEX"))
    bw.write("<!doctype html>\n")
    bw.write("<html>\n")
    makeHead(bw, TITLE)
    bw.write("<body>")
    bw.write("<h1>$TITLE</h1>")
    bw.write("<ul>")
    val dir = File("scripts")
    for (file in dir.listFiles { _, name -> name.endsWith(".sage", true) }) {
        bw.write("<li>")
        val filename = parseFile(file, INDEX, "doc")
        bw.write("<a href='$filename.html'>$filename</a>")
        bw.write("</li>")
    }
    bw.write("</ul>")
    bw.write("</body>")
    bw.write("</html>\n")
    bw.close()
}

fun makeHead(bw: BufferedWriter, title: String) {
    bw.write("<head>\n")
    bw.write("<title>$title</title>")
    bw.write("<link rel='stylesheet' type='text/css' href='style.css'>")
    bw.write("</head>\n")
}

fun parseFile(file: File, parent: String, destinationDir: String): String {
    val filename = file.name.takeWhile { it != '.' }
    val bw = BufferedWriter(FileWriter("$destinationDir/$filename.html"))
    bw.write("<!doctype html>\n")
    bw.write("<html>\n")
    makeHead(bw, filename)
    bw.write("<body>\n")
    bw.write("<a href='$parent'>Página principal</a>")
    bw.write("<h1>$filename</h1>\n")
    val br = BufferedReader(FileReader(file))
    var docType = 0
    var inClass = false
    val sb = StringBuilder()
    for (l in br.lines()) {
        if (l.isBlank()) {
            if (docType == FILE_DOC_TYPE) {
                docType = NO_DOC_TYPE
                sb.append("</p>\n")
                bw.write(sb.toString())
                sb.clear()
            }
            continue
        }
        val line = l.trim()
        if (line.startsWith("#")) {
            if (docType != NO_DOC_TYPE) {
                // Line is just #: insert line break
                if (line.length == 1) {
                    sb.append("<br>\n")
                } else {
                    val noPrefix = line.removePrefix("#").trimStart()
                    if (line.endsWith("\\\\")) {
                        sb.append("${noPrefix.removeSuffix("\\\\")}<br>\n")
                    } else {
                        sb.append(noPrefix)
                    }
                }
            } else {
                // Might be new doc
                sb.append("<p>\n")
                docType = when {
                    FILE_DOC in line -> FILE_DOC_TYPE
                    FUNCTION_DOC in line -> {
                        if (inClass) {
                            // We were inside a class and it ended
                            inClass = false
                            bw.write("</div>\n")
                            bw.write("</div>\n")
                        }
                        FUNCTION_DOC_TYPE
                    }
                    CLASS_DOC in line -> {
                        inClass = true
                        CLASS_DOC_TYPE
                    }
                    INIT_DOC in line -> INIT_DOC_TYPE
                    METHOD_DOC in line -> METHOD_DOC_TYPE
                    else -> {
                        sb.clear()
                        NO_DOC_TYPE
                    }
                }
            }
        } else if (docType != NO_DOC_TYPE) {
            // Doc ended
            sb.append("</p>\n")
            when (docType) {
                FUNCTION_DOC_TYPE, METHOD_DOC_TYPE -> {
                    if (line.startsWith("def")) {
                        val funName = line.removePrefix("def")
                                .trimStart()
                                .takeWhile { it != '(' }
                        bw.write("<div class='function'>\n")
                        bw.write("<h3>Función: $funName</h3>\n")
                        bw.write(sb.toString())
                        bw.write("</div>\n")
                    }
                }
                INIT_DOC_TYPE -> {
                    if (line.startsWith("def")) {
                        bw.write("<div class='function'>\n")
                        bw.write("<h3>Constructor</h3>\n")
                        bw.write(sb.toString())
                        bw.write("</div>\n")
                    }

                }
                CLASS_DOC_TYPE -> {
                    if (line.startsWith("class")) {
                        val className = line.removePrefix("class")
                                .trimStart()
                                .takeWhile { it != '(' && it != ':' }
                        bw.write("<div class='class'>\n")
                        bw.write("<h2>Clase: $className</h2>\n")
                        bw.write(sb.toString())
                        bw.write("<div class='methods'>\n")
                    }
                }
            }
            docType = NO_DOC_TYPE
            sb.clear()
        }
    }
    if (inClass) {
        // We were inside a class at the end of the file
        bw.write("</div>\n")
        bw.write("</div>\n")
    }
    bw.write("</body>\n")
    bw.write("</html>\n")
    br.close()
    bw.close()
    return filename
}
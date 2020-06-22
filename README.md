# Introducción al Método de los Elementos Finitos

En este curso introduciremos los conceptos básicos para resolver ecuaciones diferenciales en derivadas parciales utilizando el Método de los Elementos Finitos con herramientas de [software libre](https://es.wikipedia.org/wiki/Software_libre), dirigido a los alumnos de cursos ofrecidos por el **Grupo de Materiales Granulares** (GMG) del Departamento de Ingeniería Mecánica de la UTN - FRLP.
La idea es presentar diferentes conceptos utilizando [Jupyter notebooks](https://Jupyter.org/), y en forma activa y dinámica por lo que este curso se irá modificando según se desarrolle y a demanda de los estudiantes.

### Requerimientos

- Tener Python 3.5 o más nuevo instalado. La versión de Python puede verificarse con `python3 --version` en la línea de comandos. La última versión de Python puede descargarse de [aquí](https://www.python.org/downloads/).
- Tener [instalado Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install.html).
    - Nota: En sistemas derivados de Debian la instalación es simplemente `apt install jupyter`


### Utilización
- Clonar o descargar este repositorio.
- Ejecutar `jupyter notebook` en la línea de comando dentro del directorio del repositorio.
- Se abrirá una sesión de notebook de Jupyter en el navegador y se puede comenzar a navegar a través de los notebooks disponibles.

---

## Temas
1. [Polinomios constantes a trozo (1D)](https://nbviewer.jupyter.org/github/rirastorza/Intro2FEM/blob/master/Polinomios_constantes_atrozo/polinomios.ipynb)
2. [Elementos Finitos en 1D](https://nbviewer.jupyter.org/github/rirastorza/Intro2FEM/blob/master/Elementos_finitos_en_1D/fem1D.ipynb)<br>
    a. [FEniCS en 1D](https://github.com/rirastorza/Intro2FEM/blob/master/Elementos_finitos_en_1D/fem1D_introFEniCS.ipynb)<br>
    b. [Ejemplo mecánico](https://github.com/rirastorza/Intro2FEM/blob/master/Elementos_finitos_en_1D/mecanica1D.ipynb)<br> 
    c. [Elementos Finitos en 1D con condición de contorno Robin](https://github.com/rirastorza/Intro2FEM/blob/master/Elementos_finitos_en_1D/fem1D_Robin.ipynb)<br>
    d. Ejemplo térmico<br>
3. Elementos Finitos en 2D para problemas estacionarios, estáticos
4. Problemas dependientes del tiempo
5. Resolución iterativa
6. Mecánica de sólidos
7. Ecuación de calor
8. Mecánica de fluidos


---

## Recursos
Este curso es una recopilación de diferentes fuentes reseñadas a continuación:

1. Mats G. Larson, Fredrik Bengzon, [The Finite Element Method: Theory, Implementation, and Applications](https://www.springer.com/gp/book/9783642332869). Este libro es la referencia central del curso. Muy recomendable leerlo de principio a fin.
2. Hans Petter Langtangen, Anders Logg, [Tutorial de FEniCS](https://fenicsproject.org/tutorial/). Este es el software que utilizaremos a lo largo del curso y en el que se basarán la mayoría de los ejemplos del curso. 

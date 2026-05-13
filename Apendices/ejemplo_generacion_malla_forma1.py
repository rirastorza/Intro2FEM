#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gmsh
import os

# -------------------------
# Parámetros geométricos
L = 3.5e-3
H = 10.5e-3

# Resoluciones de malla
# (valor de Mesh.CharacteristicLengthMax)
resoluciones = {
    "gruesa": 2e-3,
    "media": 1e-3,
    "fina": 5e-4
}

# -------------------------
# Función para crear la geometría y mallar
def generar_malla(nombre_carpeta, lc_max):

    gmsh.initialize()
    gmsh.model.add("placa_3d")

    # -------------------------
    # Geometría
    # Se usa lc arbitrario en los puntos
    # porque luego se controla globalmente
    # con CharacteristicLengthMax
    lc = 1.0

    p1 = gmsh.model.geo.addPoint(-L,  L, 0, lc)
    p2 = gmsh.model.geo.addPoint( L,  L, 0, lc)
    p3 = gmsh.model.geo.addPoint( L, -L, 0, lc)
    p4 = gmsh.model.geo.addPoint(-L, -L, 0, lc)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Superficie
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Extrusión
    ext = gmsh.model.geo.extrude([(2, surf)], 0, 0, H)

    top_surface = ext[0][1]
    volume = ext[1][1]

    lateral_surfaces = [e[1] for e in ext[2:]]

    gmsh.model.geo.synchronize()

    # -------------------------
    # Physical groups

    gmsh.model.addPhysicalGroup(3, [volume], name="interno")

    gmsh.model.addPhysicalGroup(2, [top_surface], name="frente")
    gmsh.model.addPhysicalGroup(2, [surf], name="atras")

    for i, s in enumerate(lateral_surfaces):
        gmsh.model.addPhysicalGroup(2, [s], name=f"lado_{i+1}")

    # -------------------------
    # Control global de tamaño de malla

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc_max)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_max)

    # Opcional:
    # evita que Gmsh use tamaños derivados
    # de puntos o curvatura
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)

    # -------------------------
    # Generar malla
    gmsh.model.mesh.generate(3)

    # -------------------------
    # Crear carpeta
    os.makedirs(nombre_carpeta, exist_ok=True)

    # Guardar archivo
    archivo_salida = os.path.join(nombre_carpeta, "placa.msh")
    gmsh.write(archivo_salida)

    print(f"Malla guardada en: {archivo_salida}")

    gmsh.finalize()


# -------------------------
# Generar las tres mallas
for nombre, lc in resoluciones.items():
    generar_malla(nombre, lc)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

# -------------------------
# Resoluciones de malla
resoluciones = {
    "gruesa": 2e-3,
    "media": 1e-3,
    "fina": 5e-4
}

# -------------------------
# Template del .geo
geo_template = r'''
//Ejemplo 1: Placa 3d
L = 3.5e-3;
H = 10.5e-3;
Point(1) = {-L, L, 0, 1};
Point(2) = {L, L, 0, 1};
Point(3) = {L, -L, 0, 1};
Point(4) = {-L, -L, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};

Extrude {0, 0, H} {
 Surface{6};
}

Physical Volume("interno") = {1};

Physical Surface("frente",1) = {28};
Physical Surface("atras",2) = {6};

Physical Surface("abajo",3) = {27};
Physical Surface("izquierda",4) = {15};
Physical Surface("arriba",5) = {19};
Physical Surface("derecha",6) = {23};

Mesh.CharacteristicLengthMax = sizemesh;
'''

# -------------------------
# Generar las tres mallas
for nombre, sizemesh in resoluciones.items():

    # Crear carpeta
    os.makedirs(nombre, exist_ok=True)

    # Reemplazar sizemesh
    geo_content = geo_template.replace(
        "sizemesh",
        f"{sizemesh}"
    )

    # Guardar archivo .geo
    geo_file = os.path.join(nombre, "placa.geo")

    with open(geo_file, "w") as f:
        f.write(geo_content)

    # Archivo de salida
    msh_file = os.path.join(nombre, "placa.msh")

    # -------------------------
    # Ejecutar gmsh por línea de comando

    comando = f"gmsh {geo_file} -3 -o {msh_file}"

    print("Ejecutando:")
    print(comando)

    os.system(comando)

    print(f"Malla generada en: {msh_file}")

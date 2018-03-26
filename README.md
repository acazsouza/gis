
A imagem 319567_2331703_2016-12-07_0c0b-20161207T151953Z.tif é uma imagem de 
satélite multiespectral georreferenciada em formato GeoTIFF obtida pelo microsssatélite ID 0c0b da constelação PlanetScope em 7 de dezembro de 2016, às 15h19m53s UTC.

Ela possui as seguintes bandas:

ÍNDICE--------------BANDA--------------ALCANCE ESPECTRAL (nm)

===================================================

- 1-------------------Blue--------------------455 - 515
- 2-------------------Green------------------500 - 590
- 3-------------------Red---------------------590 - 670
- 4-------------------NIR---------------------780 - 860

Sua missão é calcular o seguinte:
- Percentual de área desta imagem que está coberto por algum tipo de vegetação
- Centróide geográfico da cena
- Área em quilômetros quadrados da cena
- Horário local da captura	

Esse cálculo deverá ser fornecido através de um endpoint HTTP em formato JSON.

O arquivo swagger_api.yml contém a especificação exata do formato de retorno
e do nome do endpoint.

O servidor que serve o endpoint deve ser feito em Python/Flask.

===================================================
Extra:

NDVI:
![NDVI](https://image.ibb.co/gq4Djn/ndvi.png)

MIN-X / MIN-Y:
![min-x/min-y](https://image.ibb.co/nsyg4n/minxminy.png)

MAX-X / MAX-Y:
![max-x/max-y](https://image.ibb.co/dcrjPn/maxxmaxy.png)

CENTROIDE:
![centroide](https://image.ibb.co/cpPSc7/centroide.png)
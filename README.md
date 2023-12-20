# Minimização de Erro em Ondas Parciais

O projeto tem como objetivo encontrar o ponto de erro mínimo dado uma aproximação por expansão em ondas parciais de um feixe de Bessel ideal não difrativo utilizando Beam Shape Coefficients, BSCs.

Este projeto foi desenvolvido na disciplina [SCC0713](https://gitlab.com/simoesusp/disciplinas/-/tree/master/SSC0713-Sistemas-Evolutivos-Aplicados-a-Robotica) e possui apenas fins educacionais. 

Link da apresentação: [clique aqui]().

## Motivação

Aplicação de um Algoritmo Genético a fim de comparar os resultados obtidos pela expansão de feixes em ondas parciais através de coeficientes de forma. A procura por um ponto de erro mínimo auxilia no entendimento de posicionamento de partículas submetidas, neste caso, à pressão acústica gerada pelo feixe.

## Implemetação

O algoritmo foi desenvolvido utilizando a linguagem `Julia`. O pacote `Evolutionary` para otimização evolutiva e a função `optimize` também foram utilizadas durante a implementação.

```julia
result = Evolutionary.optimize(x -> erro(x...), zeros(3), GA(populationSize = 100, selection = susinv, crossover = DC, mutation = PLM()))
```

Os parâmetros e funções são definidas abaixo:

- `Evolutionary.optimize`: função utilizada para realizar otimização evolutiva.

- `x -> erro(x...)`: função objetivo que está sendo otimizada. O operador `->` é usado para criar uma função anônima que recebe um vetor `x` e calcula `erro(x...)`.

- `zeros(3)`: ponto inicial para a otimização. Ele cria um vetor de zeros com três elementos.

- `GA(populationSize = 100, selection = susinv, crossover = DC, mutation = PLM())`:
   - `populationSize = 100`: Tamanho da população, ou seja, o número de indivíduos na população.
   - `selection = susinv`: Método de seleção, neste caso, é usado o método `susinv`.
   - `crossover = DC`: Operador de crossover, aqui é usado o operador de cruzamento discreto (`DC`).
   - `mutation = PLM()`: Operador de mutação, neste caso, é usado o operador de mutação polinomial (`PLM`).

- `result`: resultado da otimização, que inclui informações sobre a solução ótima encontrada, o valor da função objetivo nesse ponto e outras estatísticas relacionadas à execução do algoritmo evolutivo.

## Resultados

O resultado obtido foi satisfatório e auxilia a construção de um padrão estabelecido em um ponto de erro mínimo.

## Instalação

### Requisitos

Para reproduzir o programa é necessário ter instalado:

- Julia v-1.9.4
- Evolutionary
- LegendrePolynomials
- SpecialFunctions

### Execução

Para executar, basta executar os comandos abaixo na pasta do projeto:

```
julia
include("desloc_implem.jl")
using Evolutionary
result = Evolutionary.optimize(x -> erro(x...), zeros(3), GA(populationSize = 10
0, selection = susinv, crossover = DC, mutation =PLM()))
```

## Autores

- [@anajuliatagliassachi](https://github.com/anajuliatagliassachi)
- [@fda-tome](https://github.com/fda-tome)
- [@gabisandanieli](https://github.com/gabisandanieli)

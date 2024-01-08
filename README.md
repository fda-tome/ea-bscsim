# Minimização de Erro em Ondas Parciais

O projeto tem como objetivo encontrar o ponto de erro mínimo dado uma aproximação por expansão em ondas parciais de um feixe de Bessel ideal não difrativo utilizando Beam Shape Coefficients, BSCs.

Este projeto foi desenvolvido na disciplina [SCC0713](https://gitlab.com/simoesusp/disciplinas/-/tree/master/SSC0713-Sistemas-Evolutivos-Aplicados-a-Robotica) e possui apenas fins educacionais. 

Link da apresentação: [clique aqui](https://drive.google.com/file/d/1XdQnPh1N7HDIjWjLiFxOmRQNIZPdUe9t/view?usp=sharing).

## Motivação

Aplicação de um Algoritmo Genético a fim de comparar os resultados obtidos pela expansão de feixes em ondas parciais através de coeficientes de forma. A procura por um ponto de erro mínimo auxilia no entendimento de posicionamento de partículas submetidas, neste caso, à pressão acústica gerada pelo feixe.

## Implemetação

O algoritmo foi desenvolvido utilizando a linguagem `Julia`.

Os parâmetros e funções são definidas abaixo:

- `init_population(npop, bounds)`: função responsável pela inicialização da população inicial. Realizada de forma aleatória dentro dos limites do problema. 

- `crossover(parent1, parent2)`: função que realiza o cruzamento entre dois indivíduos (pais) para gerar um novo indivíduo (filho). Essa função irá calcular as médias das coordenadas x e y dos pais para gerar as coordenadas do filho.

- `mutate(individual, indpb, bounds)`: função de mutação dos indivíduos, em que é realizada uma mutação aleatória de uma porcentagem dos indivíduos de acord com o parâmetro `indpb` definido no código. A mutação realizada é pequena e dentro dos limites.

- `evaluate(population)`: função que realiza o cálculo de erro para cada indivíduo da população.

- `select(population, fitness, k, elitism_percentage)`: função de seleção dos invíduos da geração seguinte. A partir de uma porcentagem denominada `elitism_percentage` a mpopulação com o melhor fitness é copiada para a próxima geração. O restante da população é submetida ao **torneio de dois** em que são selecionados os pais da geração seguinte.

- `genetic_algorithm(error_function, ngen, npop, cxpb, indpb, bounds, elitism_percentage)`: função que gera um algoritmo genético com base nas funções anteriores e nos paramêtros utilizados.

## Resultados

O resultado obtido foi satisfatório e auxilia a construção de um padrão estabelecido em um ponto de erro mínimo.

![](/plot.png)

## Instalação

### Requisitos

Para reproduzir o programa é necessário ter instalado:

- Julia v-1.9.4
- LegendrePolynomials
- SpecialFunctions

### Execução

Para executar, basta executar os comandos abaixo na pasta do projeto:

```
julia-1.9.4 desloc_implem.jl
```

## Autores

- [@anajuliatagliassachi](https://github.com/anajuliatagliassachi)
- [@fda-tome](https://github.com/fda-tome)
- [@gabisandanieli](https://github.com/gabisandanieli)

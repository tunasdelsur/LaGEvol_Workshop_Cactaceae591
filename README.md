# Processamento de dados de sequenciamento gerados por bibliotecas de target-enrichment sequencing (Cactaceae591) - LaGEvol
# 
# 02 de junho 2023 / 04 de junho 2023
#
# Instrução: Monique, Milena e Matias
 
Esse workshop abordará algumas das principais etapas de processamento e análise de dados moleculares/genéticos provenientes de 'Targeted-enrichment Sequencing', especificamente do painel do Cactaceae591.
O curso foi organizado, de maneira práticas, com o objetivo de que você: 
  
  - 1) entenda como os dados foram e são gerados;
  - 2) aprenda a remover sequências brutas de baixa qualidade proveniente do sequenciamento;
  - 3) monte (faça o 'assemble') (d)as sequências brutas em unidades genéticas informativas (genes/locus);
  - 4) analise os arquivos de saída, identificando e removendo possíveis parálogos;
  - 5) alinhe as sequências e aprenda a ver estatísticas (e.g., % de dados faltantes - 'missing-data' -, etc.); 
  - 6) analise os alinhamentos, e filtre sequências espúrias, mal alinhadas, sem cobertura, ou ricas em gaps;
  - 7) gere árvores de máxima-verossimilhança para cada locus e para as sequências concatenadas; 
  - 8) obtenha uma árvore de espécies, utilizando um método sumário de coalescência de espécies;
 
 Além da parte prática, serão fornecidos recursos bibliográficos adicionais, e serão discutidos aspectos tangenciais à parte prática aqui aplicada.
 
 
 # 1) Filtrando e removendo sequências brutas de baixa qualidade proveniente do sequenciamento
 
Existem vários métodos e abordagens para realizar essa tarefa. Nesse tutorial, utilizaremos o _fastp_, um processador ultra-rápido e eficienete de arquivos .fastq (Veja mais aqui: https://academic.oup.com/bioinformatics/article/34/17/i884/5093234; e aqui: https://github.com/OpenGene/fastp).

*Obs.: todos os programas aqui


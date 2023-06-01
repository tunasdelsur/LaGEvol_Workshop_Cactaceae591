# Processamento de dados de sequenciamento gerados por bibliotecas de target-enrichment sequencing (Cactaceae591) - LaGEvol

![Badge em Desenvolvimento](http://img.shields.io/static/v1?label=STATUS&message=EM%20DESENVOLVIMENTO&color=GREEN&style=for-the-badge)

# 02 de junho 2023 / 04 de junho 2023

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

_*Observação_: praticamente todas as tarefas que realizaremos aqui envolvem o uso do sistema operacional Linux e da interface do _Bash_ para nos comunicar com o computador. 
O bash é um ambiente de computação que interpreta comandos de texto e executa tarefas. No sistema operacional do Windows, a partir da versão 2019, é possível instalar um subsistema do Linux no Windows (WSL).
Espera-se que todos os programas que vamos precisar para a realização desse workshop já estejam instalados e funcionando nos computadores do LaGEvol. Também é possível que os programas sejam instalados em um computador pessoal, ou em um supercomputador, mas isso não será abordado aqui.

Para interagir no bash, utiliza-se apenas linguagem de computação, que são infinitas. As principais serão brevemente exploradas e comentadas ao longo da execução de nossas tarefas. Mas, caso você nunca tenha ouvido falar nelas, ou não tenha tido nenhum contato com elas ao longo de sua vida, recomenda-se navegar na internet por outras fontes e pegar experiência com outros tutoriais: e.g., https://www.hostinger.com.br/tutoriais/comandos-linux 

## _fastp_
Para cada amostra enviada para sequenciamento, a empresa retorna dois arquivos (R1 e R2), pois é realizado um sequenciamento ao-par de cada fragmento do DNA amostrado (pair-ended).
Os arquivos vem com o nome e código do projeto de sequenciamento da empresa, com dados associados ao poço da placa em que a amostra foi montada. 
Por ex., abaixo:

```
FSC_143302_P001_WC12_132_S32_L001_R1_001.fastq.gz 

FSC_143302_P001_WC12_132_S32_L001_R2_001.fastq.gz 
```

O nome do arquivo contém:
FSC_143302 = código do projeto na empresa;
P001_WC12 = Referência ao código da placa e poço da amostra;
132_S32_L001 = Códigos dos adaptadores, e lanes usados na máquina para sequenciamento (sequenciador);
R1/R2 = Arquivo 1 ou 2 (read 1 ou read 2);
fastq.gz = código de extensão do arquivo; fastq.gz significa que a amostra está comprimida utilizando tipo de compressor gz. 

As sequencias brutas descomprimidas possuem arquivo de extensão _fastq_







A sintaxe básica para usar o fastp nos nossos





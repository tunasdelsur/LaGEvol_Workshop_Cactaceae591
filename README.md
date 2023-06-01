# Processamento de dados gerados por bibliotecas de target-enrichment sequencing (Cactaceae591)

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
 
Além da parte prática, serão fornecidos recursos bibliográficos adicionais, e discutidos aspectos tangenciais à parte prática aqui aplicada. Não abordaremos aspetos teóricos e práticos de maneira exaustiva, apenas forneceremos um panorama geral de como iniciar um projeto com dados de target-enrichment, e terminá-lo com uma árvore filogenética. Provavelmente, você irá se deparar com dificuldades teóricas, práticas e metodológicas em algum momento. Não se assuste e não desista - tod@s passam por isso, até Marie Curie, Einsten, Darwin, e Felsestein. 
 
 
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
- FSC_143302 = código do projeto na empresa;
- P001_WC12 = Referência ao código da placa e poço da amostra;
- 132_S32_L001 = Códigos dos adaptadores, e lanes usados na máquina para sequenciamento (sequenciador);
- R1/R2 = Arquivo 1 ou 2 (read 1 ou read 2);
- fastq.gz = código de extensão do arquivo; fastq.gz significa que a amostra está comprimida utilizando tipo de compressor gz. 

As sequencias brutas descomprimidas possuem arquivo de extensão _fastq_. Mas conforme avançarmos no processamento dos dados, as sequencias vão mudando de extensão, por exemplo, _.fasta_, _.phy_, _.nex_, ou outros.


A sintaxe básica para usar o _fastp_ em uma sequência é a seguinte:

```
fastp -i local/com/minha_sequencia_bruta.R1.fastq.gz -I local/com/minha_sequencia_bruta.R2.fq.gz -o local/com/sequencia_trimmada.R1.fastq.gz -O local/com/sequencia_trimmada.R2.fastq.gz
```

Com esse comando, o _fastp_ irá filtrar as sequencias de uma amostra. Porém, em um projeto de filogenômica, trabalhamos com muitos dados, e não só uma sequencia.
Então, precisamos adaptar os códigos para lidar com muitas sequências. Da mesma maneira, precisamos estar preparados para nos organizar para lidar com muitos dados. 

Você irá criar várias pastase arquivos ao longo do projeto, e você precisa estar atento para saber onde está cada uma, e qual é diferente de qual, e por quê. Ter um caderno de anotações (ou no computador, um bloco de notas - Notepad++, ou outros) pode ser útil, até estar familiarizado com o fluxo de trabalho.

Para rodar o _fastp_  em várias sequencias de uma única vez, elaboramos o seguinte código para vocês:

```
R1=(*_L001_R1_001.fastq.gz.fastq)
R2=(*_L001_R2_001.fastq.gz.fastq)
for ((i=0;i<=${#R1[@]};i++)); do fastp  -i "${R1[i]}" -I "${R2[i]}"  -o "out.${R1[i]}" -O "out.${R2[i]}" -j "${R1[i]}.fastp.json" -h "${R1[i]}.fastp.html" --dont_overwrite --failed_out "failed.${R1[i]}"; done
```

Para conferir todas as funções de um programa como o _fastp_, normalmente você pode chamá-lo no terminal com a função _--help_  ou _-h_, assim:

```
fastp -h
```

Aparecerá na sua tela todas as opções de comando do programa, com informações de cada funcionalidade.


# 2) Montando as sequências brutas (trimadas) em unidades genéticas informativas (genes/locus) com o HybPiper

## _HybPiper_

O HybPiper é um 'programa' (na verdade, chamamos de 'pipeline', pois é um programa que utiliza vários programas para chegar em um resultado final). Você pode conferir as informações completas dele aqui, inclusive de como usar: https://github.com/mossmatters/HybPiper

Ele irá juntar os as sequências brutas (raw) em "genes" ou outras unidades informativas do genoma que temos interesse. Para fazer o assembly, o HybPiper precisa dos dados brutos, um arquivo contendo as sequências das regiões alvo. além de um arquivo contendo todos os nomes das amostras, isso facilitará as próximas etapas.


1) Certifique-se de que seus dados brutos estejam desempacotados (por exemplo, sample.fastq e não sample.fastq.gz).
Se necessário, descompacte seus dados usando o comando abaixo:

```
tar -zxvf samplesdata_fastq.tar.gz
```

Verifique se você tem um arquivo contendo todos os nomes das amostras, cada uma em uma linha (ex: namelist.txt); um arquivo de referência contendo as sequências do genoma de interesse geradas com o sequenciamento do targe-capture (ex: targets.fasta); e todos os seus dados brutos trimados de amostras (ex: amostra1_R1_trimada.fastq, amostra1_R2_trimada.fastq)

Para uma amostra, usaríamos o seguinte comando:

```
hybpiper assemble -r amostra1_R1_trimada.fastq amostra1_R2_trimada.fastq -t_dna targets.fasta 
```

Novamente, lembre-se que, normalemnte, os programas tem múltiplas funções, e você pode (deve) verificar elas com o argumento _-h_:

```
hybpiper -h
```

Por exemplo, nos nossos dados, para todas as amostras, vamos executar o seguinte comando:

```
while read name; do
hybpiper assemble -r $name*.fastq -t_dna targets.fasta --prefix $name --bwa  --run_intronerate; done < namelist.txt
```



##Resumindo o enriquecimento alvo e a eficiência de recuperação de genes para um conjunto de amostras.

hybpiper stats -t_dna test_targets.fasta gene namelist.txt

#Este comando pega o resultado recuperado pelo stats e cria uma figura para visualizar a eficiência da recuperação das regiões alvo.

hybpiper recovery_heatmap seq_lengths.tsv

##Recuperando as sequências alvo.
#Esse comando busca as sequências recuperadas para o mesmo gene de várias amostras e gera um arquivo multi-FASTA desalinhado para cada gene.
 
hybpiper retrieve_sequences dna -t_dna test_targets.fasta --sample_names namelist.txt

##Recuperando sequências supercontig e/ou introns.
#Se você executou o Intronerate em suas amostras fornecendo o sinalizador --run_intronerate ao comando de montagem do hybpiper, você também pode especificar "supercontig" ou "intron" para recuperar essas sequências.

hybpiper retrieve_sequences supercontig -t_dna test_targets.fasta --sample_names namelist.txt










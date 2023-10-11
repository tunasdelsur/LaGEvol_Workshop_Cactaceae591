# Processamento de dados gerados por bibliotecas de target-enrichment sequencing (Cactaceae591)

![Badge em Desenvolvimento](http://img.shields.io/static/v1?label=STATUS&message=EM%20DESENVOLVIMENTO&color=GREEN&style=for-the-badge)

# 02 de junho 2023 / 04 de junho 2023

# Instrução: Monique, Milena, Matias
 
Esse workshop abordará algumas das principais etapas de processamento e análise de dados moleculares/genéticos provenientes de 'Targeted-enrichment Sequencing', especificamente do painel do Cactaceae591.
O curso foi organizado, de maneira teórico-prática, com o objetivo de que você: 
  
  - 1) entenda como os dados foram e são gerados;
  - 2) aprenda a remover sequências brutas de baixa qualidade proveniente do sequenciamento;
  - 3) monte (faça o 'assemble') (d)as sequências brutas em unidades genéticas informativas (genes/locus);
  - 4) analise os arquivos de saída, identificando e removendo possíveis parálogos;
  - 5) alinhe as sequências e aprenda a ver estatísticas (e.g., % de dados faltantes - 'missing-data' -, etc.); 
  - 6) analise os alinhamentos, e filtre sequências espúrias, mal alinhadas, sem cobertura, ou ricas em gaps;
  - 7) gere árvores de máxima-verossimilhança para cada locus e para as sequências concatenadas; 
  - 8) obtenha uma árvore de espécies, utilizando um método sumário de coalescência de espécies;
 
Além dessa parte, serão fornecidos recursos bibliográficos adicionais, e discutidos aspectos tangenciais à parte prática aqui aplicada. Não abordaremos aspetos teóricos e práticos de maneira exaustiva, apenas forneceremos um panorama geral de como iniciar um projeto com dados de target-enrichment, e terminá-lo com uma árvore filogenética. Provavelmente, você irá se deparar com dificuldades teóricas, práticas, metodológicas e/ou computacionais em algum momento. Não se assuste e não desista - tod@s passam por isso, até Marie Curie, Einsten, Darwin, e Felsestein. 

 
[1) Filtrando e removendo sequências brutas de baixa qualidade proveniente do sequenciamento](https://github.com/tunasdelsur/LaGEvol_course/blob/main/README.md#1-filtrando-e-removendo-sequ%C3%AAncias-brutas-de-baixa-qualidade-proveniente-do-sequenciamento)

[2) Montando as sequências brutas (trimadas) em unidades genéticas informativas (genes/locus) com o HybPiper](https://github.com/tunasdelsur/LaGEvol_course#2-montando-as-sequ%C3%AAncias-brutas-trimadas-em-unidades-gen%C3%A9ticas-informativas-geneslocus-com-o-hybpiper)

[3) Identificando e removendo possíveis parálogos](https://github.com/tunasdelsur/LaGEvol_course#3-identificando-e-removendo-poss%C3%ADveis-par%C3%A1logos)

[4) Alinhando as sequências](https://github.com/tunasdelsur/LaGEvol_course#4-alinhando-as-sequ%C3%AAncias)

[5) Filtrando sequências espúrias, mal alinhadas ou ricas em gaps e vendo as estatísticas](https://github.com/tunasdelsur/LaGEvol_Workshop_Cactaceae591/blob/main/README.md#5-filtrando-sequ%C3%AAncias-esp%C3%BArias-mal-alinhadas-ou-ricas-em-gaps-e-vendo-as-estat%C3%ADsticas)

[6) Inferência filogenética baseada em máxima-verossimilhança para cada locus e para as sequências concatenadas](https://github.com/tunasdelsur/LaGEvol_Workshop_Cactaceae591/blob/main/README.md#6-inferencia-filogen%C3%A9tica-baseada-em-m%C3%A1xima-verossimilhan%C3%A7a-para-cada-locus-e-para-as-sequ%C3%AAncias-concatenadas)

[7) Inferindo uma árvore de espécies utilizando um método sumário de coalescência de espécies](https://github.com/tunasdelsur/LaGEvol_Workshop_Cactaceae591/blob/main/README.md#7-inferindo-uma-%C3%A1rvore-de-esp%C3%A9cies-utilizando-um-m%C3%A9todo-sum%C3%A1rio-de-coalesc%C3%AAncia-de-esp%C3%A9cies)


 
 # 1) Filtrando e removendo sequências brutas de baixa qualidade proveniente do sequenciamento
 
Existem vários métodos e abordagens para realizar essa tarefa. Nesse tutorial, utilizaremos o _fastp_, um processador ultra-rápido e eficienete de arquivos `.fastq` (Veja mais aqui: https://academic.oup.com/bioinformatics/article/34/17/i884/5093234; e aqui: https://github.com/OpenGene/fastp).

_*Observação_: praticamente todas as tarefas que realizaremos aqui envolvem o uso do sistema operacional Linux e da interface do _Bash_ para nos comunicar com o computador. 
O bash é um ambiente de computação que interpreta comandos de texto e executa tarefas. No sistema operacional do Windows, a partir da versão 2019, é possível instalar um subsistema do Linux no Windows (WSL).
Espera-se que todos os programas que vamos precisar para a realização desse workshop já estejam instalados e funcionando nos computadores do LaGEvol. Também é possível que os programas sejam instalados em um computador pessoal, ou em um supercomputador, mas isso não será abordado aqui.

Para interagir no bash, utiliza-se apenas linguagem de computação, que são infinitas. As principais serão brevemente exploradas e comentadas ao longo da execução de nossas tarefas. Mas, caso você nunca tenha ouvido falar nelas, ou não tenha tido nenhum contato com elas ao longo de sua vida, recomenda-se navegar na internet por outras fontes e pegar experiência com outros tutoriais: e.g., https://www.hostinger.com.br/tutoriais/comandos-linux 

## _fastp_
Para cada amostra enviada para sequenciamento, a empresa retorna dois arquivos (`R1` e `R2`), pois é realizado um sequenciamento ao-par de cada fragmento do DNA amostrado (pair-ended).
Os arquivos vem com o nome e código do projeto de sequenciamento da empresa, com dados associados ao poço da placa em que a amostra foi montada. 
Por ex., abaixo:

```
FSC_143302_P001_WC12_132_S32_L001_R1_001.fastq.gz 

FSC_143302_P001_WC12_132_S32_L001_R2_001.fastq.gz 
```

O nome do arquivo contém:
- `FSC_143302` = código do projeto na empresa;
- `P001_WC12` = Referência ao código da placa e poço da amostra;
- `132_S32_L001` = Códigos dos adaptadores, e lanes usados na máquina para sequenciamento (sequenciador);
- `R1/R2` = Arquivo 1 ou 2 (read 1 ou read 2);
- `fastq.gz` = código de extensão do arquivo; `fastq.gz` significa que a amostra está comprimida utilizando tipo de compressor `gz`. 

As sequencias brutas descomprimidas possuem arquivo de extensão `fastq`. Mas conforme avançarmos no processamento dos dados, as sequencias vão mudando de extensão, por exemplo, `.fasta`, `.FNA`, `.phy`, `.nex`, ou outros.


A sintaxe básica para usar o _fastp_ em uma sequência é a seguinte:

```
fastp -i local/com/minha_sequencia_bruta.R1.fastq.gz -I local/com/minha_sequencia_bruta.R2.fq.gz -o local/com/sequencia_trimmada.R1.fastq.gz -O local/com/sequencia_trimmada.R2.fastq.gz
```

Com esse comando, o _fastp_ irá filtrar as sequencias de uma amostra. Porém, em um projeto de filogenômica, trabalhamos com muitos dados, e não só uma sequencia.
Então, precisamos adaptar os códigos para lidar com muitas sequências. Da mesma maneira, precisamos estar preparados para nos organizar para lidar com muitos dados. 

Você irá criar várias pastas e arquivos ao longo do projeto, e você precisa estar atento para saber onde está cada uma, e qual é diferente de qual, e por quê. Ter um caderno de anotações (ou no computador, um bloco de notas - Notepad++, ou outros) pode ser útil, até estar familiarizado com o fluxo de trabalho.

Para rodar o _fastp_  em várias sequencias de uma única vez, elaboramos o seguinte código para vocês:

```
R1=(*_L001_R1_001.fastq.gz.fastq)
R2=(*_L001_R2_001.fastq.gz.fastq)
for ((i=0;i<=${#R1[@]};i++)); do fastp  -i "${R1[i]}" -I "${R2[i]}"  -o "trim_${R1[i]}.fastq" -O "trim_${R2[i]}.fastq" -j "${R1[i]}.fastp.json" -h "${R1[i]}.fastp.html" -q 20 --dont_overwrite --failed_out "failed.${R1[i]}"; done
```

Para conferir todas as funções de um programa como o _fastp_, normalmente você pode chamá-lo no terminal com a função `--help`  ou `-h`, assim:

```
fastp -h
```

Aparecerá na sua tela todas as opções de comando do programa, com informações de cada funcionalidade.


# 2) Montando as sequências brutas (trimadas) em unidades genéticas informativas (genes/locus) com o HybPiper

## _HybPiper_

O HybPiper é um 'programa' (na verdade, chamamos de 'pipeline', pois é um programa que utiliza vários programas para chegar em um resultado final). Você pode conferir as informações completas dele aqui, inclusive de como usar: https://github.com/mossmatters/HybPiper

Ele irá juntar as sequências brutas (raw) em "genes" ou outras unidades informativas do genoma que temos interesse. Para fazer o assembly, o HybPiper precisa dos dados brutos, um arquivo contendo sequências de referência das regiões alvo, além de um arquivo contendo todos os nomes das amostras - pois isso facilitará as próximas etapas.


1) Certifique-se de que seus dados brutos estejam descompactados (por exemplo, sample.fastq e não sample.fastq.gz).
Se necessário, descompacte seus dados usando o comando abaixo:

```
tar -zxvf samplesdata_fastq.tar.gz
```

Verifique se você tem um arquivo contendo todos os nomes das amostras, cada uma em uma linha (ex: namelist.txt); um arquivo de referência contendo as sequências do genoma de interesse geradas com o sequenciamento do targe-capture (ex: targets.fasta); e todos os seus dados brutos trimados de amostras (ex: amostra1_R1_trimada.fastq, amostra1_R2_trimada.fastq)

Primeiro, ative a biblioteca do HybPiper pelo conda.

```
conda activate hybpiper
```

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
while read name; do hybpiper assemble -r $name*.fastq -t_dna targets.fasta --prefix $name --bwa  --run_intronerate; done < namelist.txt
```

Depois de fazer o assemble, vamos pegar as estatísticas, ver quantos genes foram montados para cada amostra, eficiência, etc.

```
hybpiper stats -t_dna targets.fasta gene namelist.txt
```

Agora, vamos analizar essa informação de forma visual: 

```
hybpiper recovery_heatmap seq_lengths.tsv
```

Depois de verificar a eficiência de amostra/locus, vamos sintetizar os reads recuperados para cada 'target' (gene/locus) em um alquivo multi-FASTA.

```
hybpiper retrieve_sequences dna -t_dna targets.fasta --sample_names namelist.txt --fasta_dir 01_dna
```

Além de recuperar as sequências 'target', podemos também recuperar as sequências adjacentes que foram geradas no sequenciamento. Para isso, utilizamos a o argumento "--run_intronerate" na hora de fazer o assemble, e agora recuperamos essas sequências, que o HybPiper chama de 'supercontig':

```
hybpiper retrieve_sequences supercontig -t_dna targets.fasta --sample_names namelist.txt --fasta_dir 04_supercontig
```


Ao final dessa etapa, você deverá ter como saída do HybPiper algumas pastas e arquivos, como as citadas abaixo, além de uma pasta para cada uma de suas sequências: 

```
./01_dna_seqs/
./04_supercontig_seqs/
./hybpiper_stats.tsv
./recovery_heatmap.png
./seq_lengths.tsv
```



Vamos abrir em um editor de texto um arquivo .FNA gerado pelo HybPiper e analisar o  que tem dentro dele. 


Podemos observar que além do nome da amostra, a sequência de DNA de cada amostra, temos também mensagens de aviso após o nome das amostras. Precisamos remover apenas as mensagens de erro, elas podem gerar mensagens de erro nas próximas etapas (alinhamento;trimagem) 

```
sed -i 's/*_*[0-9]*_*hits*//g' *
sed -i 's/single//g' *
sed -i 's/multi_stitched_contig_comprising_//g' *
```

# 3) Identificando e removendo possíveis parálogos

Nesse ponto, você já deve saber o que são cópias parálogas, e porque devemos estar cientes delas. Existem vários e distintos métodos de identificar e remover parálogos em conjunto de dados filogenômicas.
Nesse workshop, realizaremos um método simples, baseado no próprio HybPiper. Veja detalhes sobre isso aqui: https://github.com/mossmatters/HybPiper/wiki/Paralogs

Vamos rodar o seguinte código:

```
hybpiper paralog_retriever namelist.txt -t_dna targets.fasta
```

Depois de rodar esse código, confira o arquivo `paralog_heatmap.png`. Existem cópias parálogas no seu conjunto de dados? Quais são?
Você também pode verificar os arquivos `_paralog_report.tsv_`, `paralogs_above_threshold_report.txt`, e conferir como as cópias parálogas estão distribuídas no seu conjunto de dados e no seu grupo de estudo.

Se você tiver `mafft` e `iQTree` instalados, você pode criar uma árvore diretamente de um arquivo `*.paralogs_all.fasta` usando o seguinte comando:
```
mafft --auto gene074_paralogs_all.fasta > aligned_gene074_paralogs_all.fasta
iqtree -s aligned_gene074_paralogs_all.fasta 
```

Essas duas sequências parálogas ou alelos? 

A informação sobre cópias parálogas pode ser muito útil para identificar amostras poliplóides, duplicação de genomas, híbridos, etc. Mas para isso é necessário análises que não abordagemos no curso.

Como não veremos análises que levam em consideração cópias parálogas, é importante removê-las do nosso conjunto de dados. Dos dados do painel Cactaceae591, já temos um arquivo listando os locus que identificamos como parálogos `list_paralogs_remove.txt`.

Para remover os locos paralogos vamos primeiro criar uma pasta.
```
mkdir locos_paralogos
```

Agora podemos mover todos os arquivos fasta dos genes identificados como paralogos para a pasta locos_paralogos. Assim eles nao estaram juntos dos locos que iremos usar nas nossas análises. 

```
while read line; do mv $line ./locos_paralogos; done < list_paralogs_remove.txt
```

# 4) Alinhando as sequências
Agora que já temos as sequências dos locos e das amostras 'montadas', podemos seguir para o próximo passo, que é alinhar as sequências. Isso é necessário para que possamos comparar as sequências das amostras entre si de maneira que faça sentido evolutivo, ou seja, comparando entre sítios que tenham a mesma história evolutiva. Mesmo que para o sequênciamento de todas as amostras tenha sido usado as mesmas 'probes', é possível e provável que a sequência final montada de cada locus/amostra tenha tamanhos diferentes, tenha iniciado/terminado em regiões diferentes, ou que tenha sido montada regiões além da que fornecemos no arquivo de referência. Por isso, precisamos alinhar.

Após fazer o assemble com o HybPiper, e utilizar o comando `retrieve_sequences`, você deve ter gerado na sua pasta de trabalho os seguintes diretótios, por exemplo: `01_dna_seqs` e `04_supercontigs_seqs`.
Dentro dessas pastas, tem um arquivo `.fasta` ou `.FNA` para cada locus que foi fornecido no arquivo de referência `targest.fasta`. É cada arquivo desse que precisa ser alinhado.

## mafft

O MAFFT é um programa para alinhamento de sequências nucleotídicas ou de aminoácidos. É amplamente usado e reconhecido pela sua eficiência computacional e acurácia. Para a documentação completa, veja aqui: https://mafft.cbrc.jp/alignment/software/

O MAFFT tem diferentes algoritmos para realizar o alinhamento das sequências, que podem ser mais eficientes dependendo do tipo da sequência (longa, curta, com mais ou menos dados faltantes, etc.). Porém, ele tem uma opção para detectar 'automaticamente' qual o melhor algoritmo para alinhar as sequências quee estamos fornecendo para ele, com o argumento `--auto`.

A sintaxe básica para usar o mafft é a seguinte:

```
mafft [arguments] input > output
```

Para alinhar uma das nossas sequências montadas, poderia ser assim:

```
mafft --auto uma_sequencia.fasta > uma_sequencia_alinhada.fasta
```

Como queremos alinhar vários arquivos, na verdade todos que estão dentro da pasta, vamos utilizar um loop no bash para fazer isso.

Primeiro, vamos criar uma pasta para direcionar as sequencias alinhadas:

```
mkdir Alignments
```

Agora, o loop para alinhar todas as sequencias, e direcionar para a pasta criada:

```
nohup sh -c 'for i in *.fasta; do mafft --reorder --auto "$i" > "./Alignments/aligned_$i"; done'  &
```


# 5) Filtrando sequências espúrias, mal alinhadas ou ricas em gaps e vendo as estatísticas

Ainda que usemos os melhores programas, com os melhores algoritmos para alinhar nossos locos com múltiplas sequências, esses métodos podem falhar com certas famílias de proteínas ou em regiões específicas do alinhamento. A confiabilidade e precisão das análises dependem criticamente da qualidade dos alinhamentos, portanto se conseguirmos identificar e remover essas regiões mal alinhadas e/ou com sequências espúrias, melhoraremos a qualidade do nosso alinhamento e consequentemente as estimativas das análises que dependem da informação contida nesses alinhamentos.

Para realizar esse "polimento" dos alinhamentos utilizaremos o programa trimal. Trimal é uma ferramenta desenvolvida para fazer esse “polimento” dos alinhamentos em larga escala (genômica). De forma simplificada, ele começa lendo todas as colunas em um alinhamento e calcula uma pontuação para cada uma delas. Nesse tutorial iremos explorar dois argumentos diferentes para “polir” os alinhamentos. 
Usando o valor de gap para polir um alinhamento. 

O exemplo abaixo elimina as colunas nas quais a fração de gap for menor que 0.7 (70% das amostras).

```
trimal -in meualinhamento.fa -fasta -gt 0.7 -out meualinhamento_trim.fa -htmlout meualinhamento_trim.html 
```

Usando o argumento strict para polir o alinhamento. Esse argumento combina as informações sobre a fração de gaps em uma coluna e seus escores de similaridade.

```
trimal -in meualinhamento.fa -strict -out meualinhamento_trim.fa -htmlout meualinhamento_trim.html
```

Agora escolha um dos dois métodos e faça um loop para realizar esse polimento para todas as amostras usando um só comando. Exemplo abaixo utiliza o argumento strict:

```
for i in *.FNA; do trimal -in $i -out ./alinhamento_trimado/"$i"_trimmed.fasta -gt 0.7; done; 
```
Além do trimal, existem outros programas capazes de fazer um polimento dos dados. Um deles é o spruceup, uma ferramenta utilizada para descobrir, visualizar e remover sequências espúrias (muito discrepantes) em múltiplos alinhamentos. Essa ferramenta foi desenvolvida em Python, portanto para usa-la, além de instalar o proprio spruceup é preciso ter o Python instalado no computador. 

Para utilizar o spruceup você precisará de uma supermatrix em fasta (gerada ao concatenar os locos; para aprender como concatenar veja a etapa 6 do presente tutorial), uma árvore preliminar (opcional), e um arquivo onde todos os paramentros escolhidos para a filtragem ficarão (configuration file). Você pode encontrar um exemplo desse arquivo nos arquivos nessa página do github (my-configuration-file.conf) e editalo em qualquer editor de texto. Com todos esses arquivos prontos você já pode rodar o spruceup com o comando a seguir:

```
python -m spruceup my-configuration-file.conf
```
Após rodar essa ferramenta poderemos encontrar os seguintes arquivos de saída (output):
*0.valor-do-cut-off_nome-do-alinhamento.fasta - um arquivo fasta com seus dados trimados de acordo com os parâmetros previamente definidos (esse é seu arquivo que será usado nas próximas análises)
*report - um arquivo com a informação de quais sequências foram identificadas como espúrias e trimadas do alinhamento.
*log - um arquivo com as informaçoes que aparecem na tela enquando a análise está sendo feita.
*imagens png - São gráficos com a distribuição de distância de cada amostra. Essas figuras são interessantes de serão analisadas para confirmar que os parametros utilizados na trimagem são os mais adequados.
Para mais detalhes você pode acessar diretamente o github do autor do spruceup -> https://github.com/marekborowiec/spruceup

Agora que temos uma supermatrix com as sequencias espúrias trimadas, nós podemos separar essa supermatrix em cada loco. Lembra que ao concatenar os alinhamentos podemos pedir para o programa (que estamos utilizando para concatenar) gerar um arquivo de partição? É com esse arquivo de partição que iremos separar nossa supermatrix em locos.

```
python3 AMAS.py split -f fasta -d dna -i 0.valor-do-cut-off_nome-do-alinhamento.fasta -l partitions.txt -u fasta -j
```
Agora você verá que em sua pasta, além da supermatrix, você tem o arquivo do alinhamento 'polido' de cada um dos locos definidos no arquivo de partição.
Esses alinhamentos 'polidos" serão nossos conjunto de dados utilizados em todas as análises a seguir, portanto mova-os para uma pasta separada. 
Exemplo:
```
mkdir locos_trimados
mv OG00* ./locos_trimados
```
Pronto! Agora na sua pasta "locos_trimados" estão os locos alinhados e “polidos”.

Com esses arquivos podemos gerar algumas estatísticas para avaliar em cada loco quantas amostras temos, qual o comprimento das sequências, o número de N (caracteres indeterminados), proporção de sítios variáveis... 

Para isso utilizaremos o programa AMAS com o seguinte comando:

```
python3 AMAS.py summary -f fasta -d dna -i *.fasta -o 0.SummaryStats.csv
```

# 6) Inferencia filogenética baseada em máxima-verossimilhança para cada locus e para as sequências concatenadas
Existem diferentes programas e métodos (verossimilhança, bayesiano, coalescente...) para fazer uma inferência filogenética, aqui iremos utilizar o programa iQTree. iQTree é um programa de inferência filogenética baseado em máxima verossimilhança. De forma simplificada métodos de máxima verossimilhança vão avaliar a probabilidade de que um determinado modelo/hipótese tenha gerado os dados observados.

Para inferir uma árvore para cada loco utilizando o iQTree, podemos usar o comando abaixo:

Exemplo para um loco.
```
iqtree -s example.phy -m MFP -B 1000
```
Exemplo de loop para todos os arquivos em um só comando.
```
nohup sh -c 'for i in *.fasta; do iqtree -nt 4 -s "$i" -st DNA -m MFP -B 1000; done 2>iqtree.err' &
```

Uma outra forma de inferir uma filogenia é utilizando todos os locos de uma só vez. Para isso precisaremos gerar uma supermatrix (concatenar os dados). Concatenar os dados consiste em combinar os dados de alinhamento de sequência de vários genes em um único arquivo, arquivo esse que podemos chamar de supermatrix. 

Para concatenar os dados podemos utilizar o programa phyx, usando o comando abaixo:
```
pxcat -s *fasta -p partitions.txt -o minha_supermatrix.fasta  # -p gera um arquivo com as partições
```

Também é possível usar o amas para concatenar.
```
python3 AMAS.py concat -p partitions.txt -t minha_supermatrix.fasta -u fasta -y raxml -i *.fasta -f fasta -d dna
```

Para manter cada tipo de inferência estimada organizada, vamos criar uma pasta apenas para armazenar nossa supermatrix e as inferências geradas a partir dela. 
```
mkdir supermatrix_tree
```

Vamos mover nossa supermatrix e arquivo de partição para a pasta supermatrix_tree
```
mv minha_supermatrix.fasta ./supermatrix_tree
mv partitions.txt ./supermatrix_tree
```

Agora com a nossa supermatrix podemos utilizá-la para gerar uma árvore baseada em máxima verossimilhança utilizando todos os nossos locos.
```
iqtree -nt 4 -s minha_supermatrix.fasta  -p partitions.txt  -st DNA -m MFP -B 1000
```


# 7) Inferindo uma árvore de espécies utilizando um método sumário de coalescência de espécies.

Nessa etapa iremos inferir uma árvore de espécies. Chamamos de árvores de espécies aquelas inferências filogenéticas baseadas em coalescência. Os métodos coalescentes inferem uma filogenia incorporando a heterogeneidade genealógica esperada entre os locos (a discordância entre as árvores de genes).
 
Nesse tutorial utilizaremos o programa Astral III. Esse programa é um método sumário (sumarização de árvores) para a inferência filogenética baseado no modelo de coalescência. De maneira resumida, o programa infere a árvore de espécies que melhor concorda com a maioria dos quartetos induzidos pelas árvores de cada gene separadamente.

Como dado de entrada para o astral utilizaremos as árvores de genes geradas pelo iQTree. Para isso precisamos entrar na pasta onde estão todos as minhas árvores de genes e junta-las em um único arquivo.

```
cat *.treefile > all-gene-trees.tree
```

Pensando novamente na organização dos nossos dados e análises, vamos criar uma pasta para armazenar apenas os inputs e resultados gerados pelo Astral.
```
mkdir astral_tree
mv  all-gene-trees.tree ./astral_tree
cp ../supermatrix_tree/cnames.txt ./
cp ../supermatrix_tree/nnames.txt ./

```

Vamos renomear os tips das arvores de genes que geramos. 

```
pxrlt -t all-gene-trees.tree -c cnames.txt -n nnames.txt -o all-gene-trees.relabel.tree
```

Com esse arquivo contendo todas as nossas árvores de genes com os nós colapsados podemos rodar o programa Astral. 
Rodando o Astral: 
```
java -jar /local/da/pasta/astral.5.7.8.jar -i all-gene-trees.relabel.tree -o sptree_astral.tree 2> sptree_astral.log
```


Vamo renomear os tips da árvore de supermatrix. Voltamos então par a pasta supermatrix_tree.

```
pxrlt -t minha_supermatrix.fasta.treefile -c cnames.txt -n nnames.txt -o minha_supermatrix.fasta.relabel.tree
```

Uhuuuul!!! Agora nós já temos uma árvore de espécies para nosso conjunto de dados. :D

![Badge em Desenvolvimento](http://img.shields.io/static/v1?label=STATUS&message=EM%20DESENVOLVIMENTO&color=GREEN&style=for-the-badge)



















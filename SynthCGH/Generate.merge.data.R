# Get empirical length distribution of aCGH data (from DNAcopy segmentation of bc data)
load("length.distr.RData")
# Contiene un elenco di lunghezze più che una distribuzione, 
# Si potrebbe costruirne la cdf e fare sampling in quel modo

# Generate list structure for storing simulated data
simulated.data <- list()


# Function for generation of simulated aCGH data
generate.data <- function(le) { 
	data<-numeric()
	i<-0
	while(i<le){
 	  state<-sample(c(0,1,2,3,4,5),1,prob=c(0.01,0.08,0.81,0.07,0.02,0.01))
 	  if (state==2){  l<-sample(zero.length.distr,1)  }
 	  else { l<-sample(nonzero.length.distr,1) }
 	  if (i+l>le) { l <- le-i}
  	  data<-rbind(data,c(state,l,1))
 	  i<-i+l
	}
	return(data)
}
# Questa funzione prende in input una lunghezza (es. cromosoma)
# campiona uno stato tra 0 e 5 secondo una certa distribuzione
# e di conseguenza una lunghezza dalle liste e assembla il segmento

# In questo modo non è gestibile la simulazione
# di un pattern di alterazioni, tutti i cromosomi
# sono "uguali", nel senso che per ognuno i segmenti sono
# costruiti con le stesse distribuzioni di probabilità
# per alterazioni e conseguenti lunghezze.
# Lo stesso ragionamento potrebbe essere fatto a livello
# di bande, non è possibile forzare una certa probabilità
# su una certa banda (alterazioni segmentali specifiche,
# invece che randomiche).

# Cosa interessa a me? Poter simulare un insieme di sample
# da un Template predefinito, dove per Template si intende
# un certo tipo di alterazioni, simulando anche la 
# presenza di alterazioni casuali.
# Stato:
#	potrei avere una distribuzione di default
#	per segmenti normali e potrei specificare
#	una distribuzione manualmente per una certa
#	banda a qualsiasi livello di risoluzione
# Elimino la distribuzione delle lunghezze, 
# ma devo introdurre un possibile bias di non allineamento
# come modellarlo?
# A quel punto posso usare i vari "rumori" in questo
# generatore, più il wave effect.
# Se passo direttamente il mapping UCSC? Le lunghezze sono
# quelle.. il wave effect ha senso?
# La componente spaziale non avrebbe più senso sulla forma
# del vetrino invece che sulla vicinanza delle probe?
# Il wave effect sembrs essere dovuto al "GC content" sul genoma, 
# quindi è collegato all'ordine delle probe... non si sà perché!

# Se uso un template UCSC, conosco le probe a che intervalli
# corrispondono. Se esprimo le alterazioni in termini di
# bande, conosco su quali probe fare le estrazioni.
# Potrei farla separata per ogni probe e così simulare errori
# di ibridizzazione ed automaticamente errori nella lunghezza.
# Oppure estraggo lo stato in modo concorde, ma poi su ognuno
# mappo un errore gaussiano diverso, che può produrre lo stesso
# effetto, ma assume che il segnale sottostante sia coerente.


# Number of chromosomes
chr<-20

# Number of simulated samples
n<-500

# Generate simulated samples chr chromosomes
for (i in 1:n){
	sample<-numeric()

	# Get proportion of normal cells (copy number 2, log2 ratio = 0)
    p <- runif(1, min=0.3, max=0.7)

    # Standard deviation for noise
	sdev <- runif(1,0.1,0.2)

	for (c in 1:chr){
		matrix <- generate.data(100)
		# Per ogni cromosoma usa la funzione per assemblare il segmento
    	
    	# Get copy numbers from generated segments
		data <- as.vector(unlist(apply(matrix,1,function(x){rep(x[1],x[2])})))
		gain.loss <- rep(0,length(data))
		gain.loss[data>2] <- 1
		gain.loss[data<2] <- -1
	
    	# Get log2 ratios given a copy number change in the 1-p proportion of tumor cells
		log2ratios<-log2((data*p+2*(1-p))/2)

		# Add noise
		log2ratios<-log2ratios+rnorm(length(log2ratios),mean=0,sd=sdev)

		sample <- rbind(sample,cbind(rep(c,100),log2ratios,data,gain.loss,1:100))
	}
	sample <- cbind(sample,1:(100*chr))
	colnames(sample) <- c("Chrom","log2ratios","copynumber","gain-loss","kb","index")
		
	simulated.data[[i]] <- list(sample=sample,p=p,sdev=sdev)

}
names(simulated.data) <- paste("sample",1:n,sep="")
save(simulated.data,file="simulated.data.RData")


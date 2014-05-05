
import java.util.*;
import java.io.*;

public class Homolog{
	GenePair[] genePairs;
    GeneInfo[] genes;
	
	public Homolog(){
		genePairs = new GenePair[20];
        genes=new GeneInfo[40];
	}
    
    public Homolog(GenePair agp ){
        genePairs = new GenePair[1];
        genePairs[0] = new GenePair(agp);
        genes = new GeneInfo[2];
        genes[0] = agp.gene1;
        genes[1] = agp.gene2;
    }
    
    public Homolog[] getAllHomologs(GenomePair[] gps, int genomeNumber){
        int totalNumber = 0;
        for(int i = 0; i< gps.length; i++){
            totalNumber = totalNumber+gps[i].genePairs.length;
        }
        Homolog[] tmp = new Homolog[totalNumber];
        int index = 0;
        for(int i = 0; i< gps.length; i++){
            for(int j = 0; j< gps[i].genePairs.length; j++){
                tmp[index] = new Homolog(gps[i].genePairs[j]);
              //  tmp[index].print();
                index++;
            }
            
        }
        
        boolean moreToCombine = true;
        while(moreToCombine==true){
            moreToCombine= false;
            for(int i = 0; i< tmp.length-1; i++){
                if(tmp[i]!=null){
                    for(int j = i+1; j< tmp.length; j++){
                        if(tmp[j]!=null){
                            if(sameGroup(tmp[i],tmp[j])){
                                tmp[i] = combineTwoGroup(tmp[i],tmp[j]);
                                tmp[j] = null;
                                moreToCombine = true;
                                totalNumber--;
                            }
                            
                        }
                    }
                }
            }
        }
        Homolog[] result = new Homolog[totalNumber];
        int index1 = 0;
        for(int i = 0; i< tmp.length; i++){
            if(tmp[i]!=null){
                result[index1] =new Homolog(tmp[i]);
                result[index1].order();
                index1++;
            }
        }
        return result;
    }
    
    public boolean sameGroup(Homolog h1, Homolog h2){
		for(int i = 0; i< h1.genes.length; i++){
			for(int j = 0; j < h2.genes.length; j++){
                if(h1.genes[i].sameGene(h2.genes[j])){
                    return true;
                }
			}
		}
		return false;
	}

    
    public Homolog combineTwoGroup(Homolog h1,Homolog h2){
        Homolog result = new Homolog();
        
		GenePair[] totalgp = new GenePair[h1.genePairs.length+h2.genePairs.length];
		int indexh = 0;
		for(int i = 0; i< h1.genePairs.length; i++){
			totalgp[indexh] = h1.genePairs[i];
            indexh++;
            
		}
		for(int i = 0; i< h2.genePairs.length; i++){
			totalgp[indexh] = h2.genePairs[i];
            indexh++;
		}
		result.genePairs = new GenePair[totalgp.length];
		result.genePairs = totalgp;
        
        GeneInfo[] tmp = new GeneInfo[h1.genes.length+h2.genes.length];
        int index = 0;
        for(int i = 0; i< h1.genes.length; i++){
            tmp[index] = h1.genes[i];
            index++;
        }
        for(int i = 0; i< h2.genes.length; i++){
            if(notInList(h2.genes[i].name,tmp)){
                tmp[index] = new GeneInfo(h2.genes[i]);
                index++;
            }
        }
        result.genes = new GeneInfo[index];
        index = 0;
        for(int i = 0; i< tmp.length; i++){
            if(tmp[i]!=null){
                result.genes[index]= tmp[i];
                index++;
            }
            
        }
        return result;
		
	}
    
    public boolean notInList(String gn, GeneInfo[] list){
        for(int i = 0; i< list.length; i++){
            if(list[i]!=null && list[i].name.equals(gn)){
                return false;
            }
            if(list[i]==null){return true;}
        }
        return true;
    }
    
   	
	public void print(){
		System.out.println("a group of genes  :length ="+genePairs.length );
		/*for(int i = 0; i< genePairs.length; i++){
			if(genePairs[i]!=null){
				genePairs[i].print();
				System.out.println();
                
			}
        }*/
        for(int i = 0; i< genes.length; i++){
            genes[i].print();
        }
        System.out.println();
        System.out.println();
	}
    
    public String toStringOrthologSet(int[] allGenomeIndex){
        String[] genesForAGenome = new String[allGenomeIndex.length];
        for(int i = 0; i< genesForAGenome.length; i++){
            genesForAGenome[i] = "";
        }
        for(int i = 0; i< genes.length; i++){
            int genomeIndexhere = genes[i].genomeIndex;
            for(int j = 0; j< allGenomeIndex.length; j++){
                if(genomeIndexhere == allGenomeIndex[j]){
                    if(genesForAGenome[j].equals("")){genesForAGenome[j] =genes[i].name;}
                    else{genesForAGenome[j] =  genesForAGenome[j]+"|"+genes[i].name;}
                }
            }
        }
        for(int i = 0; i< genesForAGenome.length; i++){
            if(genesForAGenome[i].equals("")){
                genesForAGenome[i] = "missing";
            }
        }
        String result = "";
        for(int i = 0; i< genesForAGenome.length; i++){
           result = result+genesForAGenome[i]+"\t";
        }
        return result;
    }
    
    public void printAOrthologSet(int[] allGenomeIndex){
        String[] genesForAGenome = new String[allGenomeIndex.length];
        for(int i = 0; i< genesForAGenome.length; i++){
            genesForAGenome[i] = "";
        }
        for(int i = 0; i< genes.length; i++){
            int genomeIndexhere = genes[i].genomeIndex;
            for(int j = 0; j< allGenomeIndex.length; j++){
                if(genomeIndexhere == allGenomeIndex[j]){
                    if(genesForAGenome[j].equals("")){genesForAGenome[j] =genes[i].name;}
                    else{genesForAGenome[j] =  genesForAGenome[j]+"|"+genes[i].name;}
                }
            }
        }
        for(int i = 0; i< genesForAGenome.length; i++){
            if(genesForAGenome[i].equals("")){
                genesForAGenome[i] = "missing";
            }
        }
        for(int i = 0; i< genesForAGenome.length; i++){
            System.out.print(genesForAGenome[i]+"\t");
        }
            
        
    }
    
	public Homolog(Homolog ah){
        int index = 0;
        for(int i = 0; i< ah.genePairs.length; i++){
            if(ah.genePairs[i]!=null){
                index++;
            }
        }
        genePairs=new GenePair[index];
        index = 0;
		for(int i = 0; i< ah.genePairs.length; i++){
			if(ah.genePairs[i]!=null){
				genePairs[index] = new GenePair(ah.genePairs[i]);
                index++;
			}
		}
        index = 0;
        for(int i = 0; i< ah.genes.length; i++){
            if(ah.genes[i]!=null){
                index++;
            }
        }
        genes = new GeneInfo[index];
        index=0;
        for(int i = 0; i< ah.genes.length; i++){
            if(ah.genes[i]!=null){
                genes[index] = new GeneInfo(ah.genes[i]);
                index++;
            }
        }
	}
    
/*	public Homolog(GenePair[] gps){
		genePairs = new GenePair[gps.length];
        for(int i = 0; i< gps.length; i++){
            if(gps[i]!=null){
                genePairs[i] = new GenePair(gps[i]);
            }
        }
        genes = gps[0].getAllGeneList(gps);
    }*/
    
    public Homolog[] getAllHomologs1(GenePair[] gps, int genomeNumber){
        int totalNumber = gps.length;
        Homolog[] tmp = new Homolog[gps.length];
        for(int i = 0; i< gps.length; i++){
            tmp[i] = new Homolog(gps[i]);
        }
        boolean moreToCombine = true;
        while(moreToCombine==true){
            moreToCombine= false;
            for(int i = 0; i< tmp.length-1; i++){
                if(tmp[i]!=null){
                    for(int j = i+1; j< tmp.length; j++){
                        if(tmp[j]!=null){
                            if(sameGroup(tmp[i],tmp[j])){
                                tmp[i] = combineTwoGroup(tmp[i],tmp[j]);
                                tmp[j] = null;
                                moreToCombine = true;
                                totalNumber--;
                            }
                            
                        }
                    }
                }
            }
        }
        Homolog[] result = new Homolog[totalNumber];
        int index = 0;
        for(int i = 0; i< tmp.length; i++){
            if(tmp[i]!=null){
                result[index] = tmp[i];
                index++;
            }
        }
        return result;
    }
    



	
    //***********************
    //**********************
    //*******************
    public void order(){
        orderGene();
        orderGenePair();
    }
    public void orderGene(){
        for(int i = 0; i < genes.length-1; i++){
            if(genes[i]!=null){
                for(int j = i+1; j < genes.length; j++){
                    if(genes[j]!=null){
                        if(genes[i].genomeIndex>genes[j].genomeIndex){
                            switchGene(i, j);
                        }
                        if( genes[i].genomeIndex==genes[j].genomeIndex && genes[i].chrNumber>genes[j].chrNumber){
                            switchGene(i, j);
                        }
                        if(genes[i].genomeIndex==genes[j].genomeIndex && genes[i].chrNumber == genes[j].chrNumber && genes[i].positionStart > genes[j].positionStart){
                            switchGene(i, j);
                        }
                        if(genes[i].genomeIndex==genes[j].genomeIndex && genes[i].chrNumber == genes[j].chrNumber && genes[i].positionStart == genes[j].positionStart){
                            System.out.print("SPG ");
                            System.out.println( "genes[i] "+ genes[i].name + "genes[j] "+genes[j].name + " are at same chr and same startPosition");
                            // genes[j]=null;
                        }
                    }
                }
            }
		}
	}
	
	
    
	public void switchGene(int  g1, int g2){
		GeneInfo  tmp = new GeneInfo(genes[g2]);
		genes[g2] =new GeneInfo(genes[g1]);
		genes[g1] = new GeneInfo(tmp);
	}
    
		
	public void orderGenePair(){
		GenePair[] tmp = new GenePair[genePairs.length];
		for(int i = 0; i< genePairs.length; i++){
			tmp[i] = genePairs[i];
        }
		for(int i = 0; i< tmp.length-1; i++){
			for(int j = i+1; j< tmp.length; j++){
				if(tmp[i].gene1.genomeIndex>tmp[j].gene1.genomeIndex){
					GenePair agp = tmp[i];
					tmp[i] = tmp[j];
					tmp[j]=agp;
				}
				if(tmp[i].gene1.genomeIndex==tmp[j].gene1.genomeIndex  && tmp[i].gene2.genomeIndex>tmp[j].gene2.genomeIndex){
					GenePair agp = tmp[i];
					tmp[i] = tmp[j];
					tmp[j]=agp;
				}
				if(tmp[i].gene1.genomeIndex==tmp[j].gene1.genomeIndex  && tmp[i].gene2.genomeIndex==tmp[j].gene2.genomeIndex){
					int nn1 = tmp[i].gene1.positionStart;
					int nn2 = tmp[j].gene1.positionStart;
                    
                    int nn3 = tmp[i].gene2.positionStart;
					int nn4 = tmp[j].gene2.positionStart;
					if(nn1>nn2){
						GenePair agp = tmp[i];
						tmp[i] = tmp[j];
						tmp[j]=agp;
					}
                    if(nn1==nn2 && nn3>nn4){
						GenePair agp = tmp[i];
						tmp[i] = tmp[j];
						tmp[j]=agp;
					}
				}
			}
		}
		genePairs = new GenePair[tmp.length];
		genePairs = tmp;
	
	}
	
}

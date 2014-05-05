
import java.util.*;
import java.io.*;

public class GenomeInGeneInfo{
	GeneInfo[] genes;

	public GenomeInGeneInfo(GeneInfo[] g){
		genes = new GeneInfo[g.length];
		genes = g;
		int gl= 0;
		for(int i =0; i< genes.length; i++){
			if(genes[i]!=null){
				gl++;
			}
		}
	}
    
    public GenomeInGeneInfo(Homolog[] allOrthologSets, int aGenomeIndex){
        int number = 0;
        for(int i =0; i< allOrthologSets.length; i++){
            for(int j = 0; j< allOrthologSets[i].genes.length; j++){
                if(allOrthologSets[i].genes[j].genomeIndex == aGenomeIndex){
                    number++;
                }
            }
        }
        genes = new GeneInfo[number];
        int index = 0;
        for(int i =0; i< allOrthologSets.length; i++){
            for(int j = 0; j< allOrthologSets[i].genes.length; j++){
                if(allOrthologSets[i].genes[j].genomeIndex == aGenomeIndex){
                    genes[index] = new GeneInfo(allOrthologSets[i].genes[j]);
                    index++;
                }
            }
        }
        
    }
  	
	public void orderGene(){
        for(int i = 0; i < genes.length-1; i++){
            if(genes[i]!=null){
                for(int j = i+1; j < genes.length; j++){
                    if(genes[j]!=null){
                        if(genes[i]!=null && genes[j]!=null && genes[i].chrNumber>genes[j].chrNumber){
                            switchGene(i, j);
                        }
                        if(genes[i].chrNumber == genes[j].chrNumber && genes[i].positionStart > genes[j].positionStart){
                            switchGene(i, j);
                        }
                        if(genes[i].chrNumber == genes[j].chrNumber && genes[i].positionStart == genes[j].positionStart){
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
	
	
	public GenomeInString getGenomeInString(){
		int gn = 0;
		String[] tmp = new String[genes.length];
		for(int i = 0; i < tmp.length; i++){
			tmp[i] = "";
		}
		for(int i = 0; i < genes.length; i++){
			if(genes[i]!=null){
		//	here++;
				//System.out.println(
				gn++;
				int chr = genes[i].chrNumber;
				//if(chrInGenome < chr){chrInGenome = chr;}
				String aname = genes[i].newName;
				//String aname = genes[i].name;
				
				if(genes[i].orientation.equals("-")){
					aname = "-"+aname;
				}
				tmp[chr]= tmp[chr] + " " + aname;
			}
		}
		System.out.println("chr length before remove empty chr ");
		for(int i = 0; i< tmp.length; i++){
			String[] gs = splitBS(tmp[i]);
			if(gs.length>0){
				System.out.println(i+"\t"+gs.length);}
		}
		
		int chrInGenome = 0;
		for(int i = 0; i< tmp.length; i++){
			if(tmp[i].equals("")==false){
				chrInGenome++;
			}
		}
		
	//	System.out.println("         --totalChrNumber "+ chrInGenome);
	//	System.out.println("         --totalGeneNumber "+gn);
		String[] result = new String[chrInGenome];
		int indexhere = 0;
		for(int i = 0; i < tmp.length; i++){
			if(tmp[i].equals("")==false){
				result[indexhere] = tmp[i];
				//System.out.println("chr "+(indexhere+1)+"\n"+result[indexhere]);
				indexhere++;
			}
		}
        GenomeInString finalResult = new GenomeInString(result);
		return finalResult;
		
	}
	
	public String[] splitBS(String s){
		//System.out.println(s);
        String[] tmp = s.trim().split(" ");
        int index = 0;
        for(int i = 0; i < tmp.length; i++){
            if(tmp[i].trim().equals("") ==false){
                index++;
            }
        }
        String[] result = new String[index];
        index = 0;
        for(int j = 0; j <tmp.length; j++){
            if(tmp[j].trim().equals("") == false){
                result[index] = tmp[j].trim();
                index++;
            }
        }
        return result;
    }
				
}

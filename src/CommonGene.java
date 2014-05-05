// case1: first speciation. then genome doubling. order sngs, then add gray edges to try to maximize the number of same gray edges
import java.util.*;
import java.io.*;


public class CommonGene{
		
	GenomeInString[] newGenomes;
	
	public CommonGene(){
	}
	
	public void getCommonGeneForPolyploid(GenomeInString[] oldGenomes){
		//System.out.print("chrNumber  before\t");
		for(int i = 0; i< oldGenomes.length; i++){
			System.out.print(oldGenomes[i].chrs.length+"\t");
		}
		//System.out.print("geneNumber  before\t");
		for(int i = 0; i< oldGenomes.length; i++){
			int gn = countGeneNumber(oldGenomes[i].chrs);
			System.out.print(gn+"\t");
		}
		newGenomes = new GenomeInString[oldGenomes.length];
		
		String[] geneInGs = getGeneListPolyploid(oldGenomes[0].chrs);
		for(int i =1; i < oldGenomes.length; i++){
			geneInGs = updateGeneListPolyploid(oldGenomes[i].chrs, geneInGs);
		}
		for(int i = 0; i< oldGenomes.length; i++){
			String[] ng = getNewGenome(oldGenomes[i].chrs, geneInGs);
			newGenomes[i] = new GenomeInString(ng);
		}
		for(int i = 0; i< newGenomes.length; i++){
			System.out.print(newGenomes[i].chrs.length+"\t");
		}
		//System.out.print("geneNumber after\t");
		for(int i = 0; i< newGenomes.length; i++){
			int gn = countGeneNumber(newGenomes[i].chrs);
			System.out.print(gn+"\t");
		}
		
		
	}
	
	public String[] getGeneListPolyploid(String[] g){
		String[] genes = new String[80000];
		int index = 0;
		for(int i = 0; i< g.length; i++){
			String[] geneInChr = splitBS(g[i]);
			for(int j = 0; j < geneInChr.length; j++){
				String gn = geneInChr[j];
				if(gn.substring(0,1).equals("-")){
					gn = gn.substring(1);
				}
				boolean inList = false;
				for(int p = 0; p< index; p++){
					if(gn.equals(genes[p])){
						inList = true;
						break;
					}
				}
				if(inList == false){
					genes[index]= gn;
					index++;
				}
			}
		}
		return genes;
	}
	
	//geneInGs = updateGeneListPolyploid(oldGenomes[i].chrs, geneInGs);
	public String[] updateGeneListPolyploid(String[] g, String[] gl){
		String[] result = new String[gl.length];
		for(int i = 0; i< g.length; i++){
			String[] geneInChr = splitBS(g[i]);
			for(int j = 0; j < geneInChr.length; j++){
				String gn = geneInChr[j];
				if(gn.substring(0,1).equals("-")){
					gn = gn.substring(1);
				}
				int indexInList = findIndexInList(gl, gn);
				if(indexInList!=-1){
					result[indexInList] = gl[indexInList];
				}
			}
		}
		return result;
	}
	
	public int findIndexInList(String[] al, String ag){
		for(int i = 0; i< al.length; i++){
			if(al[i]!=null && ag.equals(al[i])){
				return i;
			}
		}
		return -1;
	}
	
	public void getCommonGene(GenomeInString[] oldGenomes){
		//System.out.print("chrNumber  before\t");
		for(int i = 0; i< oldGenomes.length; i++){
			System.out.print(oldGenomes[i].chrs.length+"\t");
		}
		
		//System.out.print("geneNumber  before\t");
		for(int i = 0; i< oldGenomes.length; i++){
			int gn = countGeneNumber(oldGenomes[i].chrs);
			System.out.print(gn+"\t");
		}
		newGenomes = new GenomeInString[oldGenomes.length];
		String[] geneInG1 = getGeneList(oldGenomes[0].chrs);
		int[] geneCopy = new int[geneInG1.length];
		for(int i =1; i < oldGenomes.length; i++){
			geneCopy = countCopy(geneCopy, geneInG1, oldGenomes[i].chrs); 
		}
		for(int i = 0; i< geneCopy.length; i++){
			if(geneCopy[i]<oldGenomes.length-1){
				geneInG1[i] = null;
			}
		}
		
		for(int i = 0; i< oldGenomes.length; i++){
			String[] ng = getNewGenome(oldGenomes[i].chrs, geneInG1);
			newGenomes[i] = new GenomeInString(ng);
		}
	//	System.out.print("chrNumber  after\t");
		for(int i = 0; i< newGenomes.length; i++){
			System.out.print(newGenomes[i].chrs.length+"\t");
		}
		//System.out.print("geneNumber after\t");
		for(int i = 0; i< newGenomes.length; i++){
			int gn = countGeneNumber(newGenomes[i].chrs);
			System.out.print(gn+"\t");
		}
	}
	
	public boolean hasGene(String[] gs, String ag){
		String g = ag;
		if(g.substring(0,1).equals("-")){
			g = ag.substring(1);
		}
		for(int i = 0; i< gs.length; i++){
			if(gs[i]!=null && gs[i].equals(g)){
				return true;
			}
		}
		return false;
	}
	
	public String[] removeEmptyChr(String[] ss){
		
		int chr = 0;
		for(int i = 0; i< ss.length; i++){
			if(ss[i].equals("")==false){
				chr++;
			}
		}
		String[] result = new String[chr];
		int rindex= 0;
		for(int i = 0; i< ss.length; i++){
			if(ss[i].equals("")==false){
				result[rindex] = ss[i];
				rindex++;
			}
		}
		
		return result;
	}
	
	public String[] getNewGenome(String[] og, String[] gl){
		
		String[] tmp = new String[og.length];
		for(int i = 0; i< og.length; i++){
			tmp[i] = "";
			String[] genes = splitBS(og[i]);
			for(int j = 0; j< genes.length; j++){
				if(hasGene(gl, genes[j]) == true){
				
					tmp[i] = tmp[i]+"  "+genes[j];
				}
			}
		}
		
		
		String[] result = removeEmptyChr(tmp);
		return result;
	}
	
	public int[] countCopy(int[] oldCopy, String[] geneList, String[] g){
		int[] result = new int[oldCopy.length];
		result = oldCopy;
		for(int i = 0; i < g.length; i++){
			String[] genes = splitBS(g[i]);
			for(int j = 0; j< genes.length; j++){
				String gn = genes[j];
				if(gn.substring(0,1).equals("-")==true){
					gn = gn.substring(1);
				}
				for(int k = 0; k< geneList.length; k++){
					if(gn.equals(geneList[k])){
						result[k]++;
						break;
					}
				}
			}
		}
		return result;
	}
	
	public int countGeneNumber(String[] g){
		int result = 0;
		for(int i = 0; i< g.length; i++){
			String[] genes = splitBS(g[i]);
			result = result+genes.length;
		}
		return result;
	}
	
	public String[] getGeneList(String[] g){
		String[] genes = new String[80000];
		int index = 0;
		for(int i = 0; i< g.length; i++){
			String[] geneInChr = splitBS(g[i]);
			for(int j = 0; j < geneInChr.length; j++){
				String gn = geneInChr[j];
				if(gn.substring(0,1).equals("-")){
					gn = gn.substring(1);
				}
				genes[index]= gn;
				index++;
			}
		}
		return genes;
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

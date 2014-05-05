 public class GenomePair{
	  
	  GenePair[] genePairs;
	 int genomeIndex1;
	 int genomeIndex2;
	  
	  public GenomePair(GeneInfo[] g1, GeneInfo[] g2, int[] w, int length){
		  for(int i = 0; i< g1.length; i++){
			  if(g1[i]!=null){
				  genomeIndex1 = g1[i].genomeIndex;
				  genomeIndex2 = g2[i].genomeIndex;
				  break;
			  }
		  }
          genePairs = new GenePair[length];
          for(int i = 0; i< length; i++){
              genePairs[i] = new GenePair(g1[i], g2[i],w[i]);
          }
				  
	  }
	 
  }
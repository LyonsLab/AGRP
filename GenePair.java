
public class GenePair{
	GeneInfo gene1;
	GeneInfo gene2;
	int weight;
	
	public GenePair(GeneInfo g1, GeneInfo g2, int w){
		gene1 =new GeneInfo(g1);
		gene2 = new GeneInfo(g2);
		weight = w;
		
	}
	
	
	public GenePair(GeneInfo g1, GeneInfo g2){
		gene1 =new GeneInfo(g1);
		gene2 = new GeneInfo(g2);
		weight = 0;
		
	}
	
	public GenePair(GenePair agp){
		gene1 = new GeneInfo(agp.gene1);
		gene2 = new GeneInfo(agp.gene2);
		weight = agp.weight;
		
	}
	
	public void print(){
		System.out.println("weight "+ weight);
		gene1.print();
		gene2.print();
		
	}
	
  /*  public boolean notInList(String ag, GeneInfo[] gis){
		for(int i = 0; i< gis.length; i++){
			if(gis[i]!=null){
				if(gis[i].name.equals(ag)){
					return false;
				}
			}
		}
		return true;
	}
	*/
    
    
	/*public GeneInfo[] getAllGeneList(GenePair[] gps){
		GeneInfo[] allGenes = new GeneInfo[gps.length*2];
		int index = 0;
		for(int i = 0; i< gps.length; i++){
			if(gps[i]!=null){
				String g1 = gps[i].gene1.name;
				String g2 = gps[i].gene2.name;
				if(notInList(g1, allGenes)==true){
					allGenes[index] = new GeneInfo(gps[i].gene1);
					//frequence[index] = 1;
					index++;
				}
				if(notInList(g2,allGenes)==true){
					allGenes[index] = gps[i].gene2;
					//frequence[index] = 1;
					index++;
				}
			}
		}
		GeneInfo[] genes = new GeneInfo[index];
		for(int i = 0; i< genes.length; i++){
			genes[i] = allGenes[i];
		}
		return genes;
	}*/
    
	public boolean sameGenePair(GenePair gp1, GenePair gp2){
		String gn1 = gp1.gene1.name;
		String gn2 = gp1.gene2.name;
		String gn3 = gp2.gene1.name;
		String gn4 = gp2.gene2.name;
		if(gn1.equals(gn3) && gn2.equals(gn4)){return true;}
		if(gn1.equals(gn4) && gn2.equals(gn3)){return true;}
		return false;
	}
		
}
		

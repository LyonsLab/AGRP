
public class SubgenomeRanges{
    int genomeIndex;
    int ploidy;
    int numberOfBlocks;
    int[][] ranges;
    
	public SubgenomeRanges(int gi, int p, int b, int[][] r){
        genomeIndex = gi;
        ploidy = p;
        numberOfBlocks = b;
        if(b!=0){
            ranges = new int[r.length][r[0].length];
            ranges = r;
        }else{
            ranges = new int[0][0];
        }
		
	}
    
    public SubgenomeRanges(SubgenomeRanges as){
        genomeIndex = as.genomeIndex;
        ploidy = as.ploidy;
        numberOfBlocks = as.numberOfBlocks;
        if(numberOfBlocks!=0){
        ranges = new int[as.ranges.length][as.ranges[0].length];
        for(int i = 0; i< ranges.length; i++){
            for(int j = 0; j< ranges[0].length; j++){
                ranges[i][j] = as.ranges[i][j];
            }
        }}
    }
    
    
    public SubgenomeRanges getRangesInGeneOrder(GenomeInGeneInfo ag){
        SubgenomeRanges result  = new SubgenomeRanges(this);
        for(int i = 0; i< ranges.length; i++){
            int[] p1 =getPosition(ag, this.ranges[i][2], this.ranges[i][3]);
            int[] p2 =getPosition(ag, this.ranges[i][2], this.ranges[i][4]);
            result.ranges[i][2] = p1[0];
            result.ranges[i][3] = p1[1];
            result.ranges[i][4] = p2[1];
        }
        return result;
    }
    
 /*   public int getPosition(GenomeInGeneInfo ag, int chr, int p){
        int result = -1;
        for(int i = 0; i< ag.genes.length; i++){
            if(ag.genes[i].chrNumber==chr){
                result++;
            }
            if(result==0 && p< ag.genes[i].positionStart){return 0;}
            if(i!=0 && ag.genes[i].chrNumber>chr && ag.genes[i-1].chrNumber==chr && p> ag.genes[i-1].positionEnd){return result;}
            if(ag.genes[i].chrNumber==chr && ag.genes[i].positionStart<=p && ag.genes[i].positionEnd>=p){
                return result;
                
            }
            if(i!=0 && ag.genes[i-1].chrNumber==chr && ag.genes[i].chrNumber==chr && ag.genes[i].positionStart>p && ag.genes[i-1].positionEnd<p){
                return result;
                
            }
        }
        return -1;
    }*/
    
    public int[] getPosition(GenomeInGeneInfo ag, int chr, int p){
       int[] cAndp = new int[2];
       cAndp[0] = -1; cAndp[1] = -1;
        int phere = -1;
        int resultC = 0;
        int preChere = -1;
        for(int i = 0; i< ag.genes.length; i++){
            int chere = ag.genes[i].chrNumber;  //chere 1 chr = 1
            if(chere !=preChere && chere < chr){
                resultC++;
                preChere = chere;
            }
        }
        cAndp[0] = resultC;
        int result = -1;
        for(int i = 0; i< ag.genes.length; i++){
            if(ag.genes[i].chrNumber==chr){
                result++;
            }
            if(result==0 && p< ag.genes[i].positionStart){cAndp[1] = 0; return cAndp;}
            if(i!=0 && ag.genes[i].chrNumber>chr && ag.genes[i-1].chrNumber==chr && p> ag.genes[i-1].positionEnd){cAndp[1] = result; return cAndp;}
            if(ag.genes[i].chrNumber==chr && ag.genes[i].positionStart<=p && ag.genes[i].positionEnd>=p){
                
                cAndp[1] = result; return cAndp;
                
            }
            if(i!=0 && ag.genes[i-1].chrNumber==chr && ag.genes[i].chrNumber==chr && ag.genes[i].positionStart>p && ag.genes[i-1].positionEnd<p){
                cAndp[1] = result; return cAndp;
               // return result;
                
            }
        }
        return cAndp;
    }
    
    
    public String toString(){
        String result = new Integer(genomeIndex).toString()+"\t"+ new Integer(numberOfBlocks).toString()+"\t"+new Integer(ploidy).toString()+"\n";
        result = result+"colorCode\tsubgenome\tchr\tstart\tend\n";
        for(int i = 0; i< numberOfBlocks; i++){
            for(int j = 0; j< ranges[0].length; j++){
                result = result+ranges[i][j]+"\t";
            }
             result = result+"\n";
        }
        return result;
    }
	

    
	
}
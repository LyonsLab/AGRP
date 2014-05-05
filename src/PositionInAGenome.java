  public class PositionInAGenome{
      int[][] chrs;
      int[][] positions;
      
      public PositionInAGenome(int ploidy, int geneNumber){
          chrs = new int[ploidy][geneNumber];
          positions = new int[ploidy][geneNumber];
          for(int i = 0; i< chrs.length; i++){
              for(int j = 0; j< chrs[0].length; j++){
                  chrs[i][j] = -1;
                  positions[i][j] = -1;
              }
          }
         // System.out.println("PositionForAGenome "+ chrs.length+"\t"+chrs[0].length);
      }
      
      public void addPositionForAGene(int gi, int[][] positionshere){
        //  System.out.println(gi+"\t"+positions.length+"\t"+positions[0].length);
          for(int i = 0; i<positions.length; i++){
              chrs[i][gi] = positionshere[i][0];
              positions[i][gi] = positionshere[i][1];
          }
      }
      
      
   /*   public PositionInAGenome(int ploidy, int geneNumber){
          chrs = new int[geneNumber][ploidy];
          positions = new int[geneNumber][ploidy];
          for(int i = 0; i< chrs.length; i++){
              for(int j = 0; j< chrs[0].length; j++){
                  chrs[i][j] = -1;
                  positions[i][j] = -1;
              }
          }
          // System.out.println("PositionForAGenome "+ chrs.length+"\t"+chrs[0].length);
      }
      
      public void addPositionForAGene(int gi, int[][] positionshere){
          //  System.out.println(gi+"\t"+positions.length+"\t"+positions[0].length);
          for(int i = 0; i<positions.length; i++){
              chrs[gi][i] = positionshere[i][0];
              positions[gi][i] = positionshere[i][1];
          }
      }*/
      public void print(){
          System.out.println("for a genome");
          for(int i = 0; i< chrs.length; i++){
              System.out.println("Chrcopy"+i+"\t");
              for(int j = 0; j< chrs[0].length; j++){
                  System.out.print(chrs[i][j]+"\t");
              }
              System.out.println();
              System.out.println("positioncopy"+i+"\t");
              for(int j = 0; j< chrs[0].length; j++){
                  System.out.print(positions[i][j]+"\t");
              }
              System.out.println();
          }
      }
      
}

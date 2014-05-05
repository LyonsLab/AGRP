  public class GenomeInString{
	String[] chrs;
      
	public GenomeInString(String[] s){
		chrs = new String[s.length]; 
		chrs = s;
	}
      
	  
	  public void print(){
		  for(int i = 0; i< chrs.length; i++){
			  System.out.println("chr "+i+"\n "+chrs[i]);}
	  }
	  
	  public int countGeneNumber(){
		  int result = 0;
		  for(int i = 0; i< chrs.length; i++){
			  String[] genes = splitBS(chrs[i]);
			  result = result+genes.length;
		  }
		  return result;
	  }
	  
	  public void chrLength(){
		  String[] chrLength = new String[10000];
		  for(int i = 0; i< chrLength.length; i++){
			  chrLength[i] = "";
		  }
		  for(int i = 0; i< chrs.length; i++){
			   String[] genes = splitBS(chrs[i]);
			  chrLength[genes.length]= chrLength[genes.length]+ "\t"+ new Integer(i).toString();
		  }
		  for(int i = 0; i< chrLength.length; i++){
			  if(chrLength[i].equals("")==false){
				  System.out.println("list of chrs with "+ i+"  genes"+ "\t"+chrLength[i]);
			  }
		  }
	  }
	  public String[] splitBS(String s){
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

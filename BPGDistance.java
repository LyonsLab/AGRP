import java.util.*;
import java.io.*;
public class BPGDistance{
	String[] genome1;
	String[] genome2;
	BPGPath[] paths1;
	BPGPath[] paths2;
        
        int geneNumber;
        int distance;
        
	int[] nodeInt;
	String[] nodeStr;
		String[] nodeStr1;
	String[] nodeStr2;

	/*public BPGDistance(PGMPath[] p1; PGMPath[] p2){
		//if p1 is for a leave node
		paths1 = new BPGPath[p1.length];
		//else{
		paths1 = new BPGPath[p1.
	}*/
        public BPGDistance(String[] gn1, String[] gn2){
            int index1 = 0;
            int index2 = 0;
            int gnumber1 = 0;
            int gnumber2 = 0;
            String[] tmp1 = new String[gn1.length];
            String[] tmp2 = new String[gn2.length];
            for(int i = 0; i< gn1.length; i++){
                String[] genes = splitBS(gn1[i]);
                if(genes.length!=0){
                gnumber1 = gnumber1+genes.length;
                    tmp1[index1] = gn1[i];
                    index1++;
                }
            }
		//	System.out.println(gn2.length);
              for(int i = 0; i< gn2.length; i++){
				//  System.out.println(gn2.length);

				//  System.out.println(i+"\t"+gn2[i]);
                String[] genes = splitBS(gn2[i]);
                if(genes.length!=0){
                gnumber2 = gnumber2+genes.length;
                    tmp2[index2] = gn2[i];
                    index2++;
                }
            }
          //  System.out.println(gnumber1);
			//System.out.println(gnumber2);
            if(gnumber1==gnumber2){geneNumber = gnumber1;}else{System.out.println("different gene numbers in two genomes");}
		genome1 = new String[index1];
		for(int i = 0; i < index1; i++){
                    genome1[i] = tmp1[i];
                }
                    
		genome2 = new String[index2];
		for(int i = 0; i < index2; i++){
                    genome2[i] = tmp2[i];
                }
			//System.out.println("geneNumber = "+geneNumber+ "  chr1: "+genome1.length+ "  chr2: "+genome2.length);
                
                nodeInt = new int[geneNumber*2];
                nodeStr = new String[geneNumber*2];
				nodeStr1 = new String[geneNumber*2];
			nodeStr2 = new String[geneNumber*2];
        }
	
        public void initialValue(){
            int index1 = 0;
            for(int i = 0; i < genome1.length; i++){  //w1,w2,w3,w4: weights for cs1. w6..: weights for cs2
			String[] genes = splitBS(genome1[i]);
			for(int j = 0; j < genes.length; j++){
				String fc = genes[j].substring(0,1);
				String n1 = null;
				String n2 = null;
				if(fc.equals("-")){
					n1 = genes[j].substring(1)+"h";
					n2 =  genes[j].substring(1)+"t";
					nodeInt[index1] = index1+1;
					nodeStr[index1] = n2;
				
					
					index1++;
					nodeInt[index1] =index1+1 ;
					nodeStr[index1] = n1;
					
					
					index1++;
				}else{
					n1 = genes[j]+"t";
					n2 = genes[j]+"h";
					nodeInt[index1] = index1+1;
					nodeStr[index1] = n1;
									index1++;
					nodeInt[index1] =index1+1 ;
					nodeStr[index1] = n2;
				
					
					index1++;
				}
				//System.out.print
				
			}
			}
			
			for(int i = 0; i< nodeStr1.length; i++){
				nodeStr1[i] = nodeStr[i];
				nodeStr2[i] = nodeStr[i];
			}
		//	for(int i = 0; i<nodeStr.length; i++){
		//	System.out.println( nodeInt[i] + "    "+nodeStr[i]);
		//	}
			
		//	System.out.println();
                
       //          Calendar time1 = Calendar.getInstance();
       //     System.out.println(" finish node and nodestr " + time1.getTime());	
                
			paths1 = getPath(genome1,1);
			paths2 = getPath(genome2,2);
                
            /*    for(int i = 0; i< paths1.length; i++){
                    if(paths1[i]!=null){
                        System.out.println(i + "  "+ paths1[i].toString());
                    }
                }
                
                 for(int i = 0; i< paths2.length; i++){
                    if(paths2[i]!=null){
                        System.out.println(i + "  "+ paths2[i].toString());
                    }
                }*/

                
        }
        
        public BPGPath[] getPath(String[] t1, int g){
            int nullnode = -g;
            BPGPath[] path1 = new BPGPath[geneNumber*2+1];
            for(int i = 0; i < t1.length; i++){
			String[] genes = splitBS(t1[i]);
			//if(genes.length==1){System.out.println("onegee ");}
			//g1number = g1number+genes.length;
			int preNode = 0;
			for(int j = 0; j < genes.length; j++){
				String fc = genes[j].substring(0,1);
				String n1 = null;
				String n2 = null;
				if(fc.equals("-")){
					n1 = genes[j].substring(1)+"h";
					n2 =  genes[j].substring(1)+"t";
				}else{
					n1 = genes[j]+"t";
					n2 =  genes[j]+"h";
				}
				int n1int = findNodeInt(n1,g);
				int n2int = findNodeInt(n2,g);
				if(j ==0){
					path1[n1int] = new BPGPath( n1int,nullnode,g,g);
					preNode = n2int;
					nullnode= nullnode-2;
					//t1number++;
				}
				if(j!=0 && j!=genes.length-1){
					path1[n1int] = new BPGPath(n1int, preNode,g,g);
                                        path1[preNode] = new BPGPath(preNode, n1int,g,g);
					preNode = n2int;
					//t2number++;
				}
				if(j==genes.length-1){
                                    if(genes.length!=1){
					path1[preNode] = new BPGPath(preNode, n1int,g,g);
					path1[n1int] = new BPGPath( n1int,preNode,g,g);
                                    }
					path1[n2int] = new BPGPath(n2int, nullnode,g,g);
					nullnode= nullnode-2;
				}
			}
		}
		return path1;

	}
        
	public int findNodeInt(String as, int g){
		if(g == 1){
		for(int i = 0; i < nodeStr.length; i++){
			if(nodeStr[i]!=null && nodeStr[i].equals(as)){
				nodeStr[i]=null;
				return (i+1);
			}
		}}
		if(g==2){
		for(int i = 0; i < nodeStr1.length; i++){
			if(nodeStr1[i]!=null && nodeStr1[i].equals(as)){
				nodeStr1[i]=null;
				return (i+1);
			}
		}}
		
		return 0;
	}
        
        public BPGPath getAPath(int sp, BPGPath[] ps){
            for(int i = sp; i < ps.length; i++){
                if(ps[i]!=null){
                    return new BPGPath(ps[i]);
                }
            }
            return null;
        }
	
	
		public void getValue(){
            initialValue();
         //   Calendar time = Calendar.getInstance();
         //   System.out.println(" finish initial " + time.getTime());	
            int cycleNumber=0;
            int goodPathNumber=0;
            BPGPath ap1 = getAPath(0, paths1);
            
            while(ap1!=null){
                int n1 = ap1.hn; int n2 = ap1.tn;
                int startPoint  = n1;
               // System.out.print("  "+startPoint);
                if(n1>0){paths1[n1] = null;} if(n2>0){paths1[n2] = null;}
                int nbig = n1; 
                int nsmall = n2;
                if(nbig<nsmall){nbig = n2; nsmall = n1;}
                boolean more = true;
                while( more==true){
                    if(ap1.hn<0 && ap1.tn<0){
                        more=false;   
                        if(goodCycle(ap1.hn,ap1.tn)==true){goodPathNumber++;}
                        break;
                    }
                    BPGPath l = paths2[nbig];
                    int m1 = l.hn; int m2 = l.tn;
                    if(m2>0 && m2==nsmall){ cycleNumber++; more=false; paths2[m1] = null; paths2[m2] = null;}
                    if(m2>0 && m2!=nsmall){ 
                  //   System.out.println("ap1 "+ ap1);
                        BPGPath ap2 = paths1[m2]; //ap1 = ap1.connect(ap1,ap2,l);  
                        int anotherNode = ap2.tn;
                        ap1 = new BPGPath(anotherNode, nsmall, 1,1 );
    
                        //   System.out.println("ap2 "+ ap2);
                           //   System.out.println("l "+ l);
                            //     System.out.println("newap1 "+ ap1);
                          
                        n1 = ap1.hn;  n2 = ap1.tn; nbig = n1; nsmall = n2;
                        if(nbig<nsmall){nbig = n2; nsmall = n1;}
                        int otherpa1 = paths1[m2].tn; if(otherpa1>0){paths1[otherpa1] = null;}
                        paths2[m1] = null; paths2[m2] = null;
                        paths1[m2] = null; 
                    }
                    if(m2<0 && nsmall < 0){
                        if(goodCycle(m2,nsmall)==true){
                        goodPathNumber++;}
                        more= false;
                        paths2[m1] = null;
                    }
                     if(m2<0 && nsmall > 0){
                       ap1= new BPGPath(nsmall, m2, 1, 2);
                        n1 = ap1.hn;  n2 = ap1.tn; nbig = n1; nsmall = n2;
                        if(nbig<nsmall){nbig = n2; nsmall = n1;}
                        paths2[m1] = null;
                    }
                }
                ap1 = getAPath(startPoint,paths1);
            }
            
            distance = geneNumber+genome1.length - cycleNumber-goodPathNumber;
	}
        
        public boolean goodCycle(int e1, int e2){
            if(e1>=0 || e2>=0){System.out.println("not a complete path "); return false;}
            int halfe1 = (int)e1/2;
            if(halfe1*2!= e1){ return true;}
            int halfe2 = (int)e2/2;
            if(halfe2*2!= e2){ return true;}
            return false;
            
            
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
				
					
				

   
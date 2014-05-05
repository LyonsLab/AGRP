public class BPGPath{

    int hn;
    int tn;
    int genomeForhn;
    int genomeFortn;
	   
    public BPGPath(int h, int t, int ah, int at){
        hn = h;
        tn = t;
        genomeForhn = ah;
        genomeFortn = at;        
    }
    
   public BPGPath(BPGPath ap){
		hn = ap.hn;
		tn = ap.tn;
		genomeForhn = ap.genomeForhn;
		genomeFortn = ap.genomeFortn;
    }
    
    
    public BPGPath connect(BPGPath p1, BPGPath p2, BPGPath l){  //simple version
		int h = 0; int t = 0; int gh=0; int gt = 0;
		
		if(p1.hn==l.hn && p1.genomeForhn==l.genomeForhn && p2.hn==l.tn && p2.genomeForhn==l.genomeFortn){
			h=p1.tn; gh = p1.genomeFortn; t = p2.tn; gt = p2.genomeFortn;
			if(h==0 && gh == 0){ h=p1.hn; gh = 11- p1.genomeForhn; }
			if(t==0 && gt == 0){ t = p2.hn; gt = 11-p2.genomeForhn;}
			return (new BPGPath(h,t,gh,gt));
		}
		
		if( p1.hn == l.tn  && p1.genomeForhn==l.genomeFortn && p2.hn == l.hn  && p2.genomeForhn==l.genomeForhn  ){
			h=p1.tn;  gh = p1.genomeFortn; t = p2.tn; gt = p2.genomeFortn;
			if(h==0 && gh == 0){ h=p1.hn;  gh = 11-p1.genomeForhn; }
			if(t==0 && gt == 0){  t = p2.hn; gt = 11-p2.genomeForhn; }
			return (new BPGPath(h,t,gh,gt));
		}
		
		if(p1.hn == l.tn && p2.tn == l.hn  && p1.genomeForhn==l.genomeFortn && p2.genomeFortn==l.genomeForhn){
			h = p1.tn; t = p2.hn; gh = p1.genomeFortn; gt = p2.genomeForhn;
			if(h==0 && gh == 0){ h=p1.hn;  gh =11- p1.genomeForhn; }
			if(t==0 && gt == 0){   t = p2.tn; gt =11- p2.genomeFortn;}
			return (new BPGPath(h,t,gh,gt));
		}
		
		if(  p1.hn == l.hn && p2.tn == l.tn  && p1.genomeForhn==l.genomeForhn && p2.genomeFortn==l.genomeFortn){
			h = p1.tn; t = p2.hn; gh = p1.genomeFortn; gt = p2.genomeForhn;
			if(h==0 && gh == 0){ h=p1.hn;  gh = 11-p1.genomeForhn; }
			if(t==0 && gt == 0){ t = p2.tn; gt = 11-p2.genomeFortn; }
			return (new BPGPath(h,t,gh,gt));
		}
                
		if(p1.tn==l.hn && p2.hn==l.tn  && p1.genomeFortn==l.genomeForhn && p2.genomeForhn==l.genomeFortn){ 
			h=p1.hn; t = p2.tn;gh = p1.genomeForhn; gt = p2.genomeFortn;
			if(h==0 && gh == 0){ h=p1.tn;  gh = 11-p1.genomeFortn; }
			if(t==0 && gt == 0){ t = p2.hn; gt = 11-p2.genomeForhn; }
			return (new BPGPath(h,t,gh,gt));
		}
                
		 if( p1.tn == l.tn && p2.hn == l.hn  && p1.genomeFortn==l.genomeFortn && p2.genomeForhn==l.genomeForhn){
			h=p1.hn; t = p2.tn; gh = p1.genomeForhn; gt = p2.genomeFortn;
			if(h==0 && gh == 0){ h=p1.tn;  gh = 11-p1.genomeFortn; }
			if(t==0 && gt == 0){ t = p2.hn; gt =11- p2.genomeForhn; }
			return (new BPGPath(h,t,gh,gt));
		}

		if(p1.tn==l.hn && p2.tn==l.tn  && p1.genomeFortn==l.genomeForhn && p2.genomeFortn==l.genomeFortn){
			h=p1.hn; t = p2.hn; gh = p1.genomeForhn;gt = p2.genomeForhn;
			if(h==0 && gh == 0){ h=p1.tn;  gh = 11-p1.genomeFortn; }
			if(t==0 && gt == 0){  t = p2.tn; gt =11- p2.genomeFortn; }
			return (new BPGPath(h,t,gh,gt));
		}
                
		 if( p1.tn == l.tn && p2.tn == l.hn  && p1.genomeFortn==l.genomeFortn && p2.genomeFortn==l.genomeForhn){
			h=p1.hn; t = p2.hn; gh = p1.genomeForhn;gt = p2.genomeForhn;
			if(h==0 && gh == 0){h=p1.tn;  gh = 11-p1.genomeFortn;  }
			if(t==0 && gt == 0){ t = p2.tn; gt =11- p2.genomeFortn;  }
			return (new BPGPath(h,t,gh,gt));
		}
		return null;

		
	}

    
    public String toString(){
    	String s = "a path: ";
    	if(this!=null){
        	s = "headNode is :" + hn + "("+genomeForhn+")  tailNode is "+ tn +"("+genomeFortn+")";
       	}
        return s;
    }
    
}
        
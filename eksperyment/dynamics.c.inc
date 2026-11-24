/* Read-only sampling after each physical update, independent of animation.
   Compact input cells are selected on the native, valid, non-source grid. */
static int dyn_ready=0, dyn_points=0, dyn_count[32], *dyn_cells[32];
static double dyn_interval=5, dyn_threshold=1, dyn_next=0, dyn_last=-1;
if (!dyn_ready) {
    dyn_ready=1;
    const char *spec=getenv("NEPAL_DYNAMICS_SPEC");
    if (spec) {
        FILE *input=fopen(spec,"r");
        if (!input || fscanf(input,"%d %lf %lf",&dyn_points,&dyn_interval,&dyn_threshold)!=3 ||
            dyn_points<1 || dyn_points>32 || dyn_interval<=0 || dyn_threshold<=0) exit(73);
        int total=sico.M*sico.N;
        int *lookup=malloc((size_t)total*sizeof(int));
        if (!lookup) exit(73);
        for(int k=0;k<total;k++) lookup[k]=-1;
        for(int k=0;k<sico.IMAX;k++) if(cdomain2[k]==0) lookup[px[k]*sico.N+py[k]]=k;
        for(int p=0;p<dyn_points;p++) {
            if(fscanf(input,"%d",&dyn_count[p])!=1 || dyn_count[p]<1 || dyn_count[p]>4096) exit(73);
            dyn_cells[p]=malloc((size_t)dyn_count[p]*sizeof(int));
            if(!dyn_cells[p]) exit(73);
            for(int k=0;k<dyn_count[p];k++) {
                int index;
                if(fscanf(input,"%d",&index)!=1 || index<0 || index>=total || lookup[index]<0) exit(73);
                dyn_cells[p][k]=lookup[index];
            }
        }
        free(lookup);
        if(fclose(input)) exit(73);
    }
}
if (dyn_points && tsum>dyn_last && (tsum>=dyn_next || !ccontinue || tsum>=tmax)) {
    const char *output=getenv("NEPAL_DYNAMICS_CSV");
    FILE *record=output ? fopen(output,"a") : NULL;
    if(!record || sico.MODEL!=7) exit(73);
    for(int p=0;p<dyn_points;p++) {
        int wet=0;
        double depth_sum=0, numerator=0;
        for(int k=0;k<dyn_count[p];k++) {
            int cell=dyn_cells[p][k];
            if(awt[cell][16]<dyn_threshold) continue;
            wet++;
            depth_sum+=awt[cell][16];
            numerator+=awt[cell][0]*awt[cell][13]+awt[cell][3]*awt[cell][14]+awt[cell][6]*awt[cell][15];
        }
        if(fprintf(record,"%.9f,%d,%d,%.12g,%.12g\n",tsum,p,wet,depth_sum,numerator)<0) exit(73);
    }
    if(fclose(record)) exit(73);
    dyn_last=tsum;
    dyn_next=(floor(tsum/dyn_interval)+1)*dyn_interval;
}

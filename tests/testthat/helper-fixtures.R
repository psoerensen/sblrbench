bench_fixture <- function() {
  markers<-c("m1","m2","m3");samples<-c("s1","s2","s3","s4");traits<-c("t1","t2")
  Z<-matrix(c(0,1,2,1, 1,0,1,2, 2,1,0,1),4,3,dimnames=list(samples,markers));B<-matrix(c(1,0,0,0,2,0),3,2,dimnames=list(markers,traits));G<-Z%*%B;E<-matrix(seq_len(8)/10,4,2,dimnames=list(samples,traits));Y<-G+E
  x<-list(schema_version=1L,data=list(marker_ids=markers,sample_ids=samples,trait_names=traits,train_ids=NULL,test_ids=NULL,reference_ids=NULL,genotypes=Z),truth=list(effects=B,genetic_values=G,phenotypes=Y,residuals=E,causal=list(shared="m1",specific=list(t1=character(),t2="m2"),all=c("m1","m2")),parameters=list()),scenario=list(study="test",architecture="unit",replicate=0L),provenance=list(seed=1L,simulator="fixture",transformations=character()),extras=list());class(x)<-c("sblrbench_simulation","list");x
}
oracle_result <- function(sim=bench_fixture()) new_sblrbench_result("oracle",effects=sim$truth$effects,pip=sblrbench:::.causal_matrix(sim),genetic_value=sim$truth$genetic_values)

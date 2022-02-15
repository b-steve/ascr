
#this is the function we use to interpolate mask with values from df_cov
#p is the power of weight depend on distance, when it set to be 0, it means
#all points within the search area or nearest are equally weighted;
#nn2.control is a list with controlling arguments for nn2(),
#it contains "k", "treetype", "searchtype", "radius", "eps"
#this is a lower level function, will be called by other interface function
interpolate_fun = function(mask, df_cov, p, nn2_control){
    if(is.null(nn2_control)){
        nn2_control = vector('list', 2)
        names(nn2_control) = c('data', 'query')
    }
    nn2_control$data = df_cov[,c('x', 'y')]
    nn2_control$query = mask
    nn2_output = do.call('nn2', nn2_control)


}
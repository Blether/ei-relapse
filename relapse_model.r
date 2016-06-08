# model relapse data and predict longer term outcomes
# data from EI discharge 2013 audit with relapse count
# 2016-04

nozeroes_poissonmodel = function(realdata=data) {
	# 	ignore zero count data as excess expected and run quick poisson model to get starting estimate of prob
	data = realdata
	datalength = length(data$n_relapses)
	
	model = glm(count[2:datalength] ~ n_relapses[2:datalength], data = data, family='poisson')
	return (model)
}

NBfit_NBsize_zeroignore = function(size, prob, counts=data$count, relapses=data$n_relapses) {
	# use with optimize() to find best fit value of size parameter given prob
	# doesn't use the zero relapses value
	dist_from_model = dnbinom(seq(0, max(relapses)), size, prob)
	total_people = sum(counts)
	total_model_zi = sum(dist_from_model[2:length(dist_from_model)])
	counts_from_model = dist_from_model[2:length(dist_from_model)]*total_people/total_model_zi
	differences = counts_from_model - counts
	resids = sum(differences^2)
	return (resids)
}

NBfit_NBprob_zeroignore = function(prob, size, counts=data$count, relapses=data$n_relapses) {
	# use with optimize() to find best fit value of prob parameter given size
	# doesn't use the zero relapses value
	dist_from_model = dnbinom(seq(0, max(relapses)), size, prob)
	total_people = sum(counts)
	total_model_zi = sum(dist_from_model[2:length(dist_from_model)])
	counts_from_model = dist_from_model[2:length(dist_from_model)]*total_people/total_model_zi
	differences = counts_from_model - counts
	resids = sum(differences^2)
	return (resids)
}

NBfit_optimise_zeroignore = function(NB_size, NB_prob, counts, relapses, size_bounds = c(.5,2), prob_bounds=c(.3,1)){
	# 	alternate fitting of size and prob for negative binomial model until stable value of size emerges
	# doesn't use the zero relapses value
	size_change = 1
	prob_change = 1
	i = 0
	while ( abs(size_change) > .001 )  {
		old_NB_size = NB_size
		old_NB_prob = NB_prob
		NB_size = optimise(NBfit_NBsize_zeroignore, size_bounds, prob = NB_prob, counts[2:length(counts)], relapses[2:length(relapses)])$minimum
		NB_prob = optimise(NBfit_NBprob_zeroignore, prob_bounds, size = NB_size, counts[2:length(counts)], relapses[2:length(counts)])$minimum
		size_change = (NB_size - old_NB_size)/NB_size
		prob_change = (NB_prob - old_NB_prob)/NB_prob
		i = i + 1
		}
	return (list(NB_prob = NB_prob, NB_size = NB_size, size_change = size_change, prob_change = prob_change, iterations_required = i))
}

relapse_model_fit = function(data, NB_prob='', NB_size='') {
	# 	fits a negative binomial model to relapse data and estimates the excess number without any relapse
	
	if ( NB_prob == '' ) {
	# for initial model fit ignore the zero relapse group as has excess zeroes
	crudepmodel = nozeroes_poissonmodel(data)
	NB_prob = exp(coef(crudepmodel)[2])[[1]]
	NB_size = 1 # size initial value set arbitrarily at present
	
	NB_new_fit = NBfit_optimise_zeroignore(NB_size, NB_prob,
	 counts=data$count, relapses=data$n_relapses)
	NB_prob = NB_new_fit$NB_prob
	NB_size = NB_new_fit$NB_size
	} 
	
	# get zero excess given best fit model
	dist_from_model = dnbinom(seq(0, max(data$n_relapses)), NB_size, NB_prob)
	total_people_zi = sum(data$count[2:length(data$count)])
	total_model_zi = sum(dist_from_model[2:length(dist_from_model)])
	scale_factor = total_people_zi/total_model_zi
	
	zero_excess_abs = data$count[1]-round(dnbinom(0, NB_size, NB_prob)*scale_factor, 0)
	zero_excess_prop = zero_excess_abs/sum(data$count)
	counts_from_model = dnbinom(data$n_relapses, NB_size, NB_prob)*(total_people_zi + (data$count[1] - zero_excess_abs))
	
	# 	for goodness of fit, don't include the zero-events value which is modelled not really observed - more conservative
	e_for_chisq = counts_from_model[2:length(counts_from_model)]
	o_for_chisq = data$count[2:length(data$count)]
# 	df_for_chisq = (length(data$count)-1)-1-1
# 	chistat = sum((o_for_chisq - e_for_chisq)^2/e_for_chisq)
# 	fit_p = pchisq(chistat, df_for_chisq, lower.tail = FALSE) # where p-value is large indicates no evidence of lack of fit
	# 	consider whether can use given small n in each bin at times - maybe ks.test() instead
	fit_p =  ks.test(cumsum(o_for_chisq), cumsum(e_for_chisq))$p
	
	output = list(size = NB_size, prob = NB_prob, counts_from_model = counts_from_model, zero_excess_prop = zero_excess_prop, zero_excess_abs = zero_excess_abs, n_relapses = data$n_relapses, rawcounts = data$count, fit_p = fit_p)
	return (output) 
}

plot_model_fit = function(outputlist, ...) {
	# plot from the output of relapse_model_fit()
	counts = outputlist$rawcount
	zero_excess = outputlist$zero_excess_prop
	n_relapses = outputlist$n_relapses
	counts_model = outputlist$counts_from_model
	NBprob = outputlist$prob
	NBsize = outputlist$size
	abs_best_estimate = outputlist$zero_excess_abs
	fit_p = outputlist$fit_p

	where = barplot(counts, names.arg=n_relapses, xlab='number of relapses', ylab = 'number of people', main = 'distribution of numbers of relapses with EI services', border=0, col='grey', ...)
	barplot(counts[1], add=1, border=0, axes=0, col='light blue')
	barplot((counts[1]-abs_best_estimate), add=1, border=0, axes=0, col='pink')
	arrows(where[1]+.6, counts[1], where[1]+.6, counts[1]-abs_best_estimate, length=.05, code=3)
	abs_of_zero_group = abs_best_estimate/counts[1]
	label_text = paste('ultra low risk (', 100*round(zero_excess, 2), '% of the total', '; ', 100*round(abs_of_zero_group,2), '% of those without relapse)', sep='')
	ze_label_y = counts[1]-(abs_best_estimate/2)
	text(where[1]+.6, ze_label_y, label_text, pos=4)

	lines(where[1:2], counts_model[1:2], type='l', col='red', lty='dashed')
	segments(where[1]-.5, counts[1]-abs_best_estimate, where[1]+.5, col='red', lty='dashed')
	lines(where[2:length(where)], counts_model[2:length(where)], type='l', col='red')	
	
	text_about_NB = paste('majority fit with negative binomial model\n', '(size ', round(NBsize, 2), ', prob ', round(NBprob, 2), '; p ', round(fit_p, 2), ')', sep='')
	text(where[5], counts[2], text_about_NB)
	return (where)
}

adjust_quantile_for_zero_excess = function(quantile, N, zero_excess_abs) {
	adjquantile = ((N*quantile)-zero_excess_abs)/(N - zero_excess_abs)
	return(adjquantile)
}

future_model_predictions = function(max_n_relapse, size, prob, N, zero_excess_abs, extrapolation_factor) {
	range_n = seq(0, max_n_relapse)
	newsize = size * extrapolation_factor
	counts_from_model = dnbinom(range_n, newsize, prob)*(N - zero_excess_abs)
	zero_relapse_prop = (zero_excess_abs + round(counts_from_model[1], 0))/N
	ten_or_more_prop = 1 - pnbinom(10, newsize, prob)
	median_predicted = qnbinom(adjust_quantile_for_zero_excess(.5, N, zero_excess_abs), newsize, prob)
	IQR_predicted = c(qnbinom(adjust_quantile_for_zero_excess(.25, N, zero_excess_abs), newsize, prob), qnbinom(adjust_quantile_for_zero_excess(.75, N, zero_excess_abs), newsize, prob))
	
	output=list(zero_excess_abs = zero_excess_abs, counts_from_model = counts_from_model, zero_relapse_prop = zero_relapse_prop, ten_or_more_prop = ten_or_more_prop, median_predicted = median_predicted, IQR_predicted = IQR_predicted, prob = prob, newsize = newsize, range_n = range_n)
	return (output)
}

summary_text_future_model = function(future_model) {
	outtext = paste('no relapses: ',
	 round(future_model$zero_relapse_prop, 2)*100, 
	 '%\nten or more relapses: ', round(future_model$ten_or_more_prop, 2)*100, 
	 '%\nmedian relapses: ', future_model$median_predicted, 
	 ' (IQR ', future_model$IQR_predicted[1], '-', future_model$IQR_predicted[2], ')',
	  sep='')
	return (outtext)
}

plot_future_model = function (future_model, ...) {
	range_n = future_model$range_n
	zero_excess = future_model$zero_excess
	counts_from_model = future_model$counts_from_model
	foo = counts_from_model
	foo[1] = foo[1] + zero_excess
	where = barplot(foo, names.arg=range_n, xlab='number of relapses', border=0, col='pink', main='model predictions for future', ...)
	barplot(foo[1], add=1, border=0, axes=0, col='light blue')
	barplot(counts_from_model[1], add=1, border=0, axes=0, col='pink')
	
	annotation = paste('At 10 years:\n', summary_text_future_model(future_model),sep='')
	text(where[8], max(foo), annotation)
	
	return (where)
}

double_analysis = function(data, duration, label, future_duration = 10) {
	outputmodel = relapse_model_fit(data)
	future_model = future_model_predictions(10, outputmodel$size, outputmodel$prob, sum(outputmodel$rawcounts), outputmodel$zero_excess_abs, (future_duration/duration))

	par(bty='n')
	par(mfcol=c(2,1))
	where = plot_model_fit(outputmodel, sub=label)
	a = plot_future_model(future_model, ylim=c(0,max(outputmodel$rawcounts)))

	return (outputmodel)
}


teamdata = read.csv('relapse_by_team.csv')

step_duration_days = 755
leo_duration_days = 804
coast_duration_days = 1067
teamtotal_duration_days = 863
leostep_duration_days = 776

leostep_output = double_analysis(data.frame(n_relapses = teamdata$n_relapses, count = (teamdata$LEO + teamdata$STEP)), leostep_duration_days/365, 'LEO+STEP', 10)

# slamoutput = relapse_model_fit(data.frame(n_relapses = teamdata$n_relapses, count = teamdata$total))

leo_output = relapse_model_fit(data.frame(n_relapses = teamdata$n_relapses, count = teamdata$LEO), NB_prob = leostep_output$prob, NB_size = leostep_output$size) 
step_output = relapse_model_fit(data.frame(n_relapses = teamdata$n_relapses, count = teamdata$STEP), NB_prob = leostep_output$prob, NB_size = leostep_output$size) 
coast_output = relapse_model_fit(data.frame(n_relapses = teamdata$n_relapses, count = teamdata$COAST), NB_prob = leostep_output$prob, NB_size = leostep_output$size ) 

par(mfcol=c(3,1))
plot_model_fit(leo_output, sub='LEO')
plot_model_fit(step_output, sub='STEP')
plot_model_fit(coast_output, sub='COAST')



# set.seed(421)

# 設定100樣本、5個變數，第一個變數為黃金標準變數。
# n = 100; p = 5;
# r.x = matrix(rnorm(n * p), , p) # 隨機取樣原始資料。
# r.z = r.x[ ,1] + rnorm(n) # 選取第一個變數為黃金變數。
# x = scale(r.x) # 將raw data標準化。
# z = scale(r.z) # 將第一個變數標準化。
# h = n^(-1 / 2)
# l = c(1, .2, .4, .3, .5)
# ind.d.l = c(.1, .05, 0, .03, 0)

# ------------------------------------------------

# testing function example:
# dscrt(x, z, l)
# cntin(x, z, l, h)
# d.theta.sh.h.p(x, z, l , h)
# optimal.delta(x, z, l, h, ind.d.l)
# (t1 = cgAUC(r.x, r.z, h, delta = 1, auto = FALSE, tau = 1)) # the delta be constant
# (t2 = cgAUC(r.x, r.z, h, delta = 1, auto = TRUE, tau = 1)) # the delta be variable

# ------------------------------------------------

# Smooth function.
# t: a value, the difference between any two subjects.
# h: a value, the function width. n^(-1/2 ~ -1/5), or =.3 in constant

s.h = function(t, h)
{
	1 / (1 + exp(-t / h))
}

# ------------------------------------------------

# original method
# discrete function, single variable
# y: a vactor(subjects)
# z: a vactor

dscrt = function(y, z, l) # theta_sh_0
{
	n = length(z) # 找出幾個subjects
	psii = rep(0, n)
	temp = 0
	y = y %*% l

	for(i in 1:n)
	{
		psij = 0
		for(j in 1:n)
		{
			if(i != j)
			{
				temp = ifelse(((y[i] - y[j]) * (z[i] - z[j])) > 0, 1, ifelse(((y[i] - y[j]) == 0 || (z[i] - z[j]) == 0), 0.5, 0))
				psij = psij + temp
			}
		}
		psii[i] = psij # 有n個值
	}
	theta.h.p = sum(psii) / (n * (n - 1)) # formula 9
	
	var = sum(((psii / (n - 1)) - theta.h.p)^2) / ((n / 2)*((n / 2) - 1)) # formula 10 & 11
	
	return(list(theta.h.p = theta.h.p, var = var))
}

# ------------------------------------------------

cntin = function(y, z, l, h) # theta_sh
{
	n = dim(y)[1]
	psii = rep(0, n)
	temp = 0

	for(i in 1:n)
	{
		psij = 0
		for(j in 1:n)
		{
			if(i != j)
			{
				temp = s.h(t(y[i, ] - y[j, ]) %*% l, h) * (2 * s.h((z[i] - z[j]), h) - 1) + 1 - s.h((z[i] - z[j]), h)
				psij = psij + temp
			}
		}
		psii[i] = psij
	}
	theta.sh.h.p = sum(psii) / (n * (n - 1))
	var = sum(((psii / (n - 1)) - theta.sh.h.p)^2) / ((n / 2)*((n / 2) - 1))

	return(list(theta.sh.h.p = theta.sh.h.p, var = var))
}

# ------------------------------------------------

d.theta.sh.h.p = function(y, z, l, h) # 微分, step 4
{
	n = dim(y)[1]
	sij = temp = 0
	d.theta.sh.h.p = rep(0, dim(y)[2])

	for(i in 1:n)
	{
		for(j in 1:n)
		{
			if(i != j)
			{
				temp = (s.h(t(y[i, ] - y[j, ]) %*% l, h) * (1 - s.h(t(y[i, ] - y[j, ]) %*% l, h)) * (2 * s.h(z[i] - z[j], h) - 1) * (y[i, ] - y[j, ])) / h
				sij = sij + temp
			}
		}
	}
	d.theta.sh.h.p = sij / (n * (n - 1))
	return(d.theta.sh.h.p)
}

# ------------------------------------------------

optimal.delta = function(y, z, l, h, ind.d.l)
{
	l.i = matrix(rep(l, times = 50), nrow = 50, byrow = TRUE)
	
	delta = seq(0, 5, length = 50)
	m = delta %*% t(ind.d.l) # 增加的量

	l.i = l.i + m

	l.i.max = apply(l.i, 1, max)
	l.i = l.i / l.i.max

	theta = rep(0, 50) # AUC = theta
	
	for(i in 2:50)
	{
		theta[i] = cntin(y, z, l.i[i, ], h)$theta.sh.h.p
	}
	
	delta.star = delta[which(theta == max(theta))]

	return(delta.star)
}

# ------------------------------------------------

cgAUC = function(x, z, h, delta = 1, auto = FALSE, tau = 1)
{
	x = scale(x);
	z = scale(z);
	
	conv = FALSE # 是否收斂，設定為「否」
	n = dim(x)[1] # 取出X矩陣的列數目
	p = dim(x)[2] # 取出X矩陣的行數目
	cntin.ri = dscrt.ri = rep(0, p)
	id = diag(p)

	for(i in 1:p)
	{
		dscrt.ri[i] = dscrt(x, z, id[i, ])$theta.h.p
		cntin.ri[i] = cntin(x, z, id[i, ], h)$theta.sh.h.p
	}

	beta.i = ifelse(cntin.ri > 0.5, 1, -1) # 校正AUC的變數方向

	dscrt.ri = ifelse(dscrt.ri > 0.5, dscrt.ri, (1 - dscrt.ri))
	cntin.ri = ifelse(cntin.ri > 0.5, cntin.ri, (1 - cntin.ri))

	y = x * matrix(beta.i, n, p, byrow = TRUE) # 調整過後的x
	max.x = which(cntin.ri == max(cntin.ri))

	theta.sh.h.p = 0
	l = id[max.x, ]

	while(conv == FALSE)
	{
		d.l = d.theta.sh.h.p(y, z, l, h)
		max.d.l = max(d.l)
		ind.d.l = ifelse(d.l >= (tau * max.d.l), 1, 0) * d.l
		
		if (auto == TRUE)
		{
			delta = optimal.delta(y, z, l, h, ind.d.l)
		}

		l = l + delta * ind.d.l
		l = l / max(l)
		theta.temp = cntin(y, z, l, h)$theta.sh.h.p

		ifelse(abs(theta.temp - theta.sh.h.p) < 0.0001, conv <- TRUE, conv <- FALSE) # 僅能用箭頭
		theta.sh.h.p = theta.temp
	}

	optimal.dscrt = dscrt(y, z, l)
	theta.sh.h.p.var = cntin(y, z, l, h)$var

	l = l * beta.i

	return(list(
	l = l,

	theta.sh.h.p = theta.sh.h.p,
	theta.sh.h.p.var = theta.sh.h.p.var,
	cntin.ri = cntin.ri,

	theta.h.p = optimal.dscrt$theta.h.p,
	theta.h.p.var = optimal.dscrt$var,
	dscrt.ri = dscrt.ri,
	
	delta = delta
	))
}

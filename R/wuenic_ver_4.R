library(zoo)
library(rolog)

ccode = "jam"
args = commandArgs(trailingOnly=TRUE)
if(length(args))
    ccode = tools::file_path_sans_ext(args[1])

once(call("load_files", sprintf("countries/%s.pl", ccode), list(call("encoding", quote(iso_latin_1)))))

Vn = c("bcg", "bcgx", "dtp1", "dtp1x", "dtp3", "dtp3x", 
       "hepb0", "hepb1", "hepb3", "hepb3x", "hepbb","hepbbx", "hib1", "hib3", "hib3x",
       "opv1", "ipv1", "ipv1x", "ipv2", "ipv2_frac", "ipv2x", "mcv1", "mcv1x", "mcv2", "pcv1", "pcv3", "pcv3x",
       "pol1", "pol3", "pol3x", "rcv1", "rotac", "rotacx", "yfv", "menga")
Yn = 1985:2022

# MG, discuss: khm.pl has an "opv1" vaccine, I guess it is ipv1

sawtooth = 10
svy.thrs = 10
svy.scope = 2
unpd.thrs = 10

YV = list(Y=Yn, V=Vn)
YV_bool = matrix(FALSE, nrow=length(Yn), ncol=length(Vn), dimnames=YV)
YV_int = matrix(NA_integer_, nrow=length(Yn), ncol=length(Vn), dimnames=YV)
YV_real = matrix(NA_real_, nrow=length(Yn), ncol=length(Vn), dimnames=YV)
YV_char = matrix(NA_character_, nrow=length(Yn), ncol=length(Vn), dimnames=YV)

# Simulate Prolog's "outward" rounding
tround <- function(number)
  sign(number) * trunc(abs(number) + 0.5 + sqrt(.Machine$double.eps))

# Default setting
#
# firstRubellaAtSecondMCV(_C, rcv1, _Y, mcv2).

firstRubellaAtSecondMCV = YV_char
firstRubellaAtSecondMCV[, "rcv1"] = "mcv2"

# Prolog atoms to R character strings, variables to R character strings, list
# elements of type name:elem to named list elements

atom2char = function(q)
{
  if(is.expression(q))
    return(iconv(as.character(q), "latin1", "latin1"))

  if(is.symbol(q))
    return(iconv(as.character(q), "latin1", "latin1"))
  
  # Translate mmr/rubella to string
  if(is.call(q))
    if(as.character(q[[1]]) == "/")
      return(paste(collapse="/", as.character(q[-1])))

  if(is.call(q))
  { args <- as.list(q)
    args[-1] <- lapply(args[-1], FUN=atom2char)
    return(as.call(args))
  }

  if(is.list(q) && !is.null(names(q)))
    return(lapply(q, atom2char))
  
  if(is.list(q))
  { n = NULL
    for(i in 1:length(q))
    { item = as.list(q[[i]])
      q[[i]] = atom2char(item[[3]])
      n = c(n, atom2char(item[[2]]))
    }
    names(q) = n
    return(q)
  }
  
  return(q)
}

# Read all pred(C, Y, Out) and save Out in a vector
rep3 = function(pred="births_UNPD")
{ q = call(pred, expression(C), expression(Y), expression(Out))
  s = findall(q)
  s = lapply(s, atom2char)
  s = lapply(s, as.data.frame)
  s = do.call("rbind", s)
  s$Y = as.character(s$Y)

  m = NA + numeric(length(Yn))
  names(m) = Yn
  m[s$Y] = s$Out
  return(m)
}

# Read all pred(C, V, Y, Out) and save Out in a matrix
rep4 = function(ccode, pred="admin")
{
  if(!is.list(once(call("current_predicate", call("/", as.symbol(pred), 4L)))))
    return(YV_int)
  
  q = call(pred, as.name(ccode), expression(V), expression(Y), expression(Out))
  s = findall(q)
  s = lapply(s, atom2char)
  s = lapply(s, as.data.frame)
  s = do.call("rbind", s)
  s$Y = as.character(s$Y)

  index = which(!(s$V %in% Vn))
  if(length(index))
  { warning("Unknown vaccine(s): ", paste(unique(s$V[index]), collapse=", "))
    s = s[-index, ]
  }

  m = YV_int
  m[cbind(s$Y, s$V)] = s$Out
  return(m)
}

est_req = function()
{
  q = call("estimate_required", expression(C), expression(V), expression(Y), 
           expression(Comb), expression(A5))
  s = findall(q)
  s = lapply(s, atom2char)
  s = lapply(s, as.data.frame)
  s = do.call("rbind", s)
  s$Y = as.character(s$Y)
  
  index = which(!(s$V %in% Vn))
  if(length(index))
  { warning("Unknown vaccine(s): ", paste(unique(s$V[index]), collapse=", "))
    s = s[-index, ]
  }

  m = YV_bool
  m[cbind(s$Y, s$V)] = TRUE
  return(m)
}

rubella = function()
{
  q = call("estimate_required", expression(C), expression(V), expression(Y), 
           expression(Comb), expression(Rub))
  s = findall(q)
  s = lapply(s, atom2char)
  s = lapply(s, as.data.frame)
  s = do.call("rbind", s)
  s$Y = as.character(s$Y)
  
  m = YV_char
  m[cbind(s$Y, s$V)] = s$Rub # Revert back
  return(m)
}

# Multiple surveys per year, therefore no YV-matrix, but a long data frame
survey_results = function(ccode)
{
  if(!is.list(once(call("current_predicate", quote(survey_results/6L)))))
    return(data.frame(V=NULL, Y=NULL, Yn=NULL, Id=NULL, Info=NULL, Cov=NULL))
  
  q = call("survey_results", as.name(ccode), expression(V), expression(Y), 
    expression(Id), expression(Info), expression(Cov))
  s = findall(q)
  s = lapply(s, atom2char)
  s = lapply(s, as.data.frame)
  s = do.call("rbind", s)
  s$Yn = as.character(s$Y)

  index = which(!(s$V %in% Vn))
  if(length(index))
  { warning("Unknown vaccine(s): ", paste(unique(s$V[index]), collapse=", "))
    s = s[-index, ]    
  }

  return(s)
}

# Multiple decisions per year, therefore no YV-matrix, but a long data frame
wgd = function()
{
  # Convenience function: Replace NULL by NA in lists
  padNA = function(q, v)
  { qq = q[v]
    names(qq) = v
    qq[] = ifelse(lapply(qq, is.null), NA, qq)
    return(qq)
  }
  
  q = call("wgd", expression(C), expression(V), expression(Y0), expression(Y1),
           expression(Dec), expression(Info),
           expression(Id), expression(Cov), expression(A3), expression(A4))
  s = findall(q)
  s = lapply(s, atom2char)
  s = lapply(s, padNA,
    v=c("C", "V", "Y0", "Y1", "Dec", "Info", "Id", "Cov", "A3", "A4"))
  s = lapply(s, as.data.frame)
  s = do.call("rbind", s)

  # Fix strange year ranges
  s$Y0 = pmax(s$Y0, min(Yn))
  s$Y1 = pmin(s$Y1, max(Yn))
  
  # V = NA means that a decision applies to all vaccines
  V.na = function(d)
  { if(is.na(d['V']))
      return(data.frame(V=Vn, Y0=d['Y0'], Y1=d['Y1'], Dec=d['Dec'], 
        Id=d['Id'], Info=d['Info'], Cov=d['Cov'], row.names=NULL))

    data.frame(V=d['V'], Y0=d['Y0'], Y1=d['Y1'], Dec=d['Dec'],
      Id=d['Id'], Info=d['Info'], Cov=d['Cov'])
  }

  s = apply(s, MARGIN=1, simplify=FALSE, FUN=V.na)
  s = do.call("rbind", s)
  
  # If a year range is given, apply decision to each included year
  Y.range = function(d)
  { data.frame(V=d['V'], Y=d['Y0']:min(2022, d['Y1']),
      Dec=d['Dec'], Id=d['Id'], Info=d['Info'], Cov=d['Cov'], row.names=NULL)
  }

  s = apply(s, MARGIN=1, simplify=FALSE, FUN=Y.range)
  s = do.call("rbind", s)
  s$Cov = as.integer(s$Cov)
  s$Y = as.character(s$Y)
  
  # The string "NA" does not need to be shown
  # s$Info[which(s$Info == "NA")] = ""
  return(s)
}

# 1. Load country-specific information
s = once(call("country", expression(Code), expression(Country)))
Code = atom2char(s$Code)
Country = atom2char(s$Country)

s = once(call("date", expression(Date)))
Date = atom2char(s$Date)

Ereq = est_req()
Rubella = rubella()
Admin = rep4(ccode, "admin")
Gov = rep4(ccode, "gov")
Legacy = rep4(ccode, "legacy")
Vaccinated = rep4(ccode, "vaccinated")
Target = rep4(ccode, "target")
Births = rep3("births_UNPD")[as.character(Yn)]
Surviving = rep3("si_UNPD")[as.character(Yn)]
Survey = survey_results(ccode)
Decisions = wgd()

# 2. Add information from government and administration.
#
# Reported to WHO and UNICEF is government estimate. If government
# estimate missing, then reported is administrative data. If both
# missing, fail.
#
# R implementation is in reverse order, with the higher-order predicate
# eventually overwriting the lower predicate.

Rep.Cov = YV_real
Rep.Src = YV_char

# reported(C, V, Y, Source, Coverage) :-
#     admin0(C, V, Y, Cov0),
#     (   decision(C, V, Y, ignoreGov, _, _, _)
#     ;   not(gov(C, V, Y, _))
#     ),
#     not(decision(C, V, Y, ignoreAdmin, _, _, _)),
#     !,
#     Source = admin,
#     Coverage = Cov0.

Rep.Cov = Admin
Rep.Src = ifelse(is.na(Admin), NA, "admin")

# Working group decides to ignore admin data
ignore = Decisions$Dec == "ignoreAdmin"
index = cbind(Decisions$Y[ignore], Decisions$V[ignore])
Rep.Cov[index] = NA
Rep.Src[index] = NA

# Use admin data only if government data is invalid
# Todo: This code can be removed, data will be overwritten below

gov = !is.na(Gov)
ignore = (Decisions$Dec == "ignoreGov")
gov[cbind(Decisions$Y[ignore], Decisions$V[ignore])] = FALSE
Rep.Cov[gov] = NA
Rep.Src[gov] = NA

# reported(C, V, Y, Source, Coverage) :-
#    gov(C, V, Y, Cov0),
#    not(decision(C, V, Y, ignoreGov, _, _, _)),
#    !,
#    Source = gov,
#    Coverage = Cov0.

gov = !is.na(Gov)
ignore = (Decisions$Dec == "ignoreGov")
gov[cbind(Decisions$Y[ignore], Decisions$V[ignore])] = FALSE
Rep.Cov[gov] = Gov[gov]
Rep.Src[gov] = "gov"

# 3. Check reported data
#
# Reasons to exclude reported data are working group decisions, coverage > 100%
# or temporal inconsistency.

reject = YV_bool

# reported_rejected(C, V, Y) :-
#     reported(C, V, Y, _, Coverage),
#     not(reported_later(C, V, Y)),
#     Prec is Y - 1,
#     reported(C, V, Prec, _, PrecCov),
#     sawtooth_threshold(Threshold),
#     (   member(V, [pcv3, rotac])
#     ->  PrecCov - Coverage > Threshold
#     ;   abs(PrecCov - Coverage) > Threshold
#     ), !.

# Sudden decline in most recently reported data for new vaccines
V.new = Rep.Cov[, Vn %in% c("pcv3", "rotac"), drop=FALSE]

# Search for last reported data
Diff = apply(V.new, 2, na.trim, sides="right", simplify=FALSE) # last reported
Diff = lapply(Diff, diff, simplify=FALSE)                      # jumps
Diff = lapply(Diff, rev)                                       # -> first
Diff = lapply(Diff, `[`, 1)                                    # jump of interest

J = sapply(Diff, `<`, -sawtooth)                               # check for decline
Y = sapply(Diff, names)[which(J)]
V = names(Diff)[which(J)]
reject[cbind(Y, V)] = J[which(J)]

Rep.Prev = rbind(Rep.Cov, NA)
Rep.Prev[] = rbind(NA, Rep.Cov)

Rej.Info = YV_char
Rej.Info[] = ""
Rej.Info[cbind(Y, V)] = sprintf(
    "Reported data excluded due to decline in reported coverage from %i level to %i percent. ",
    Rep.Prev[cbind(Y, V)], Rep.Cov[cbind(Y, V)])

# Sudden change in most recently reported data for classic vaccines
V.new = Rep.Cov[, !(Vn %in% c("pcv3", "rotac")), drop=FALSE]
Diff = apply(V.new, 2, na.trim, sides="right", simplify=FALSE)
Diff = lapply(Diff, diff)
Diff = lapply(Diff, rev)
Diff = lapply(Diff, `[`, 1)
Diff = lapply(Diff, abs)                                     # up or down

J = sapply(Diff, `>`, sawtooth)
Y = sapply(Diff, names)[which(J)]
V = names(Diff)[which(J)]
reject[cbind(Y, V)] = J[which(J)]

Rep.Prev = rbind(Rep.Cov, NA)
Rep.Prev[] = rbind(NA, Rep.Cov)

Rej.Info[cbind(Y, V)] = sprintf(
    "Reported data excluded due to sudden change in coverage from %i level to %i percent. ",
    Rep.Prev[cbind(Y, V)], Rep.Cov[cbind(Y, V)])

# reported_rejected(C, V, Y) :-
#     reported(C, V, Y, _, Coverage),
#     Prec is Y - 1,
#     Succ is Y + 1,
#     reported(C, V, Prec, _, PrecCov),
#     reported(C, V, Succ, _, SuccCov),
#     sawtooth_threshold(Threshold),
#     (   Coverage - PrecCov > Threshold,
#         Coverage - SuccCov > Threshold
#     ;   PrecCov - Coverage > Threshold,
#         SuccCov - Coverage > Threshold
#     ), !.

up    = apply(rbind(NA, Rep.Cov), 2, diff) > sawtooth
down  = apply(rbind(Rep.Cov, NA), 2, diff) < -sawtooth
index = which(up & down, arr.ind=TRUE)
reject[index] = TRUE

prev = index
prev[, 1] = prev[, 1] - 1
succ = index
succ[, 1] = succ[, 1] + 1
Rej.Info[index] = sprintf(
    "Reported data excluded due to an increase from %i percent to %i percent with decrease to %i percent. ",
    Rep.Cov[prev], Rep.Cov[index], Rep.Cov[succ])

down  = apply(rbind(NA, Rep.Cov), 2, diff) < -sawtooth
up    = apply(rbind(Rep.Cov, NA), 2, diff) > sawtooth
index = which(down & up, arr.ind=TRUE)
reject[index] = TRUE

prev = index
prev[, 1] = prev[, 1] - 1
succ = index
succ[, 1] = succ[, 1] + 1
Rej.Info[index] = sprintf(
    "Reported data excluded due to decline in reported coverage from %i percent to %i percent with increase to %i percent. ",
    Rep.Cov[prev], Rep.Cov[index], Rep.Cov[succ])

# % Implausible coverage
# reported_rejected(C, V, Y) :-
#     reported(C, V, Y, _, Coverage),
#     Coverage > 100,
#     !.

index = which(Rep.Cov > 100, arr.ind=TRUE)
reject[index] = TRUE
Rej.Info[index] = sprintf(
    "Reported data excluded because %i percent greater than 100 percent. %s", 
    Rep.Cov[index], Rej.Info[index])

# reported_rejected(C, V, Y) :-
#     decision(C, V, Y, acceptReported, _, _, _),
#     !,
#     fail.

index = Decisions[Decisions$Dec == "acceptReported", ]
reject[cbind(index$Y, index$V)] = FALSE
Rej.Info[cbind(index$Y, index$V)] = NA

# % Ignore dominates accept
# reported_rejected(C, V, Y) :-
#     decision(C, V, Y, ignoreReported, _Expl0, _, _),
#     !.

index = Decisions[Decisions$Dec == "ignoreReported", ]
reject[cbind(index$Y, index$V)] = TRUE

# Rep.Cov[reject] = NA
# Rep.Src[reject] = NA

# 4. Time series of reported data
#
# % Reported data available
# reported_time_series(C, V, Y, Source, Coverage) :-
#     estimate_required(C, V, Y, _, _),
#     reported(C, V, Y, Source0, Cov0),
#     not(reported_rejected(C, V, Y)),
#     !,
#     Source = Source0,
#     Coverage = Cov0.

TS.Cov = YV_int
TS.Src = YV_char

TS.Cov[!reject & Ereq] = Rep.Cov[!reject & Ereq]
TS.Src[!reject & Ereq] = Rep.Src[!reject & Ereq]

# index = !Ereq
# TS.Cov[!index] = NA
# TS.Src[!index] = NA

# % Interpolation, no data/reported data excluded between two years
# reported_time_series(C, V, Y, Source, Coverage) :-
#     estimate_required(C, V, Y, _, _),
#     (   not(reported(C, V, Y, _, _))
#     ;   reported_rejected(C, V, Y)
#     ),
#     year_before_reported(C, V, Y, Prec, PrecCov),
#     year_after_reported(C, V, Y, Succ, SuccCov),
#     !,
#     Source = interpolated,
#     interpolate(Prec, PrecCov, Succ, SuccCov, Y, Coverage).

Rep = Rep.Cov
Rep[reject] = NA
inter = apply(Rep, 2, na.approx, na.rm=FALSE)
index = Ereq & (is.na(Rep.Cov) | reject) & !is.na(inter)
TS.Cov[index] = tround(inter[index])
TS.Src[index] = "interpolated"

# % Extrapolation, latest required estimate
# reported_time_series(C, V, Y, Source, Coverage) :-
#     estimate_required(C, V, Y, _, _),
#     (   not(reported(C, V, Y, _, _))
#     ;   reported_rejected(C, V, Y)
#     ),
#     nearest_reported(C, V, Y, _Year, Cov0),
#     !,
#     Source = extrapolated,
#     Coverage = Cov0.

Rep = Rep.Cov
Rep[reject] = NA
extra = apply(Rep, 2, zoo::na.locf, na.rm=FALSE)

# Don't overwrite interpolated
index = Ereq & (is.na(Rep.Cov) | reject) & is.na(inter) & !is.na(extra)
TS.Cov[index] = round(extra[index])
TS.Src[index] = "extrapolated"

Rep = Rep.Cov
Rep[reject] = NA
extra = apply(Rep, 2, zoo::na.locf, na.rm=FALSE, fromLast=TRUE)
index = Ereq & (is.na(Rep.Cov) | reject) & is.na(inter) & !is.na(extra)
TS.Cov[index] = round(extra[index])
TS.Src[index] = "extrapolated"

# 5. Survey data
#
# % Survey results passed for inclusion in the analysis include:
# % card or history results for cohorts 12-23, 18-29, 15-26, 24-35 months of age
# survey_for_analysis(C, V, Y, ID, Description, Coverage) :-
#     survey_results0(C, V, Y, ID, Description, Coverage),
#     member(confirm:'card or history', Description),
#     member(age:AgeCohort, Description),
#     member(AgeCohort, ['12-23 m', '18-29 m', '15-26 m', '24-35 m']).

Svy.Cov = YV_int
Svy.Expl = YV_char
Svy.Info = array(NA_character_, dim=c(length(Yn), length(Vn), 0),
                 dimnames=list(Y=Yn, V=Vn, Id=NULL))

cnf = Survey$Info.confirm == "card or history"
age = Survey$Info.age %in% c("12-23 m", "18-29 m", "15-26 m", "24-35 m")

Dn = levels(as.factor(Survey$Id))
if(any(cnf & age))
{
  index    = Survey[cnf & age, ]
  
  # MG, discuss: Only needed to find Ids for surveys that are ignored for multiple
  # reasons
  # Example: bcg,1999,afg2000231 ignored because of sample size as well as by a
  # wgd
  Svy.CoH = array(NA_integer_, dim=c(length(Yn), length(Vn), length(Dn)), 
                  dimnames=list(Y=Yn, V=Vn, Id=Dn))
  Svy.CoH[cbind(index$Y, index$V, index$Id)] = index$Cov
  
  # % Reasons to exclude a survey include:
  # %    Sample size < 300,
  # %    The working group decides to exclude the survey.
  # survey_accepted(C, V, Y, ID, Coverage) :-
  #     survey_for_analysis(C, V, Y, ID, Desc, Cov0),
  #     (   decision(C, V, Y, acceptSurvey, _, ID, _)
  #     ;   member(ss:Size, Desc),
  #         Size >= 300,
  #         not(decision(C, V, Y, ignoreSurvey, _, ID, _)),
  #         not(decision(C, V, Y, ignoreSurvey, _, na, _))
  #     ),
  #     % Check if survey needs to be modified
  #     (   survey_modified(C, V, Y, ID, _, Modified)
  #     ->   Coverage = Modified
  #     ;   Coverage = Cov0
  #     ).
  
  index    = Survey[cnf & age, ]
  
  Dn = levels(as.factor(Survey$Id))
  Svy.Title = array(NA_character_, dim=c(length(Yn), length(Vn), length(Dn)), 
                    dimnames=list(Yn, Vn, Dn))
  Svy.Title[cbind(index$Y, index$V, index$Id)] = index$Info.title
  
  index    = Survey[cnf & age, ]
  Svy.Ana = array(NA_integer_, dim=c(length(Yn), length(Vn), length(Dn)), 
                  dimnames=list(Y=Yn, V=Vn, Id=Dn))
  Svy.Ana[cbind(index$Y, index$V, index$Id)] = index$Cov
  
  # % Recall bias is estimated by comparing the first and third dose of a vaccine
  #
  # vaccine(dtp3, dtp1).
  # vaccine(pol3, pol1).
  # vaccine(hib3, hib1).
  # vaccine(hepb3, hepb1).
  # vaccine(pcv3, pcv1).
  #
  # survey_modified(C, V, Y, ID, Expl, Coverage) :-
  #     member(V, [dtp3, pol3, hib3, hepb3, pcv3]),
  #
  #     % Third dose, card only
  #     survey_results0(C, V, Y, ID, DescriptionCard3Dose, C3Cov),
  #     member(confirm:card, DescriptionCard3Dose),
  #     member(age:AgeCohortCard3Dose, DescriptionCard3Dose),
  #     member(AgeCohortCard3Dose, ['12-23 m', '18-29 m', '15-26 m', '24-35 m']),
  
  V13 = c(dtp1="dtp3", pol1="pol3", hib1="hib3", hepb1="hepb3", pcv1="pcv3")
  vac = Survey$V %in% V13
  cnf = Survey$Info.confirm == "card"
  age = Survey$Info.age %in% c("12-23 m", "18-29 m", "15-26 m", "24-35 m")
  index = Survey[vac & cnf & age, ]
  
  Svy.C3 = array(NA_integer_, dim=c(length(Yn), length(V13), length(Dn)), 
                 dimnames=list(Y=Yn, V=V13, Id=Dn))
  Svy.C3[cbind(index$Y, index$V, index$Id)] = index$Cov
  
  #     % First dose, card or history
  #     vaccine(V, First),
  #     survey_results0(C, First, Y, ID, DescriptionCoH1Dose, CoH1Cov),
  #     member(confirm:'card or history', DescriptionCoH1Dose),
  #     member(age:AgeCohortCoH1, DescriptionCoH1Dose),
  #     member(AgeCohortCoH1, ['12-23 m', '18-29 m', '15-26 m', '24-35 m']),
  
  vac = Survey$V %in% names(V13)
  cnf = Survey$Info.confirm == "card or history"
  age = Survey$Info.age %in% c("12-23 m", "18-29 m", "15-26 m", "24-35 m")
  index = Survey[vac & cnf & age, ]
  
  Svy.CH1 = array(NA_integer_, dim=c(length(Yn), length(names(V13)), length(Dn)), 
                  dimnames=list(Y=Yn, V=names(V13), Id=Dn))
  Svy.CH1[cbind(index$Y, index$V, index$Id)] = index$Cov
  
  #     % First dose, card only
  #     survey_results0(C, First, Y, ID, DescriptionCard1Dose, C1Cov),
  #     C1Cov > 0,
  #     member(confirm:card, DescriptionCard1Dose),
  #     member(age:AgeCohortCard1Dose, DescriptionCard1Dose),
  #     member(AgeCohortCard1Dose, ['12-23 m', '18-29 m', '15-26 m', '24-35 m']),
  
  vac = Survey$V %in% names(V13)
  cnf = Survey$Info.confirm == "card"
  age = Survey$Info.age %in% c("12-23 m", "18-29 m", "15-26 m", "24-35 m")
  cv0 = Survey$Cov > 0
  index = Survey[vac & cnf & age & cv0, ]
  
  Svy.C1 = array(NA_integer_, dim=c(length(Yn), length(names(V13)), length(Dn)), 
                 dimnames=list(Y=Yn, V=names(V13), Id=Dn))
  Svy.C1[cbind(index$Y, index$V, index$Id)] = index$Cov
  
  #     Adj is C3Cov / C1Cov,
  #     ThirdHistoryAdj is (CoH1Cov - C1Cov) * Adj,
  #     CovAdjusted is C3Cov + ThirdHistoryAdj,
  #     bound_0_100(CovAdjusted, Cov0),
  
  Adj = Svy.C3 / Svy.C1
  H3Adj = Adj * (Svy.CH1 - Svy.C1)
  H3Adj[] = pmin(99, pmax(0, tround(Svy.C3 + H3Adj)))
  
  #     survey_for_analysis(C, V, Y, ID, Description, SurveyCoverage),
  #     Cov0 \= SurveyCoverage,
  
  # MG, discuss: Should we use some minimum distance instead of rounding?
  index = which(tround(Svy.Ana[, V13, , drop=FALSE]) != H3Adj, arr.ind=TRUE)
  index1 = cbind(Y=Yn[index[, "Y"]], V=names(V13)[index[, "V"]], Id=Dn[index[, "Id"]])
  index3 = cbind(Y=Yn[index[, "Y"]], V=V13[index[, "V"]], Id=Dn[index[, "Id"]])
  
  #     SurveyCovRounded is round(SurveyCoverage),
  #     CH1 is round(CoH1Cov),
  #     C1 is round(C1Cov),
  #     C3 is round(C3Cov),
  #     member(title:Title, Description),
  #     concat_atom([Title, ' card or history results of ', SurveyCovRounded,
  #         ' percent modified for recall bias to ', Cov0,
  #         ' percent based on 1st dose card or history coverage of ',
  #         CH1, ' percent, 1st dose card only coverage of ',
  #         C1, ' percent and 3rd dose card only coverage of ',
  #         C3, ' percent. '], Expl),
  #         Coverage = Cov0.
  
  Svy.Info = array(NA_character_, dim=c(length(Yn), length(Vn), length(Dn)),
                   dimnames=list(Y=Yn, V=Vn, Id=Dn))
  
  if(any(index))
  {
    Svy.Info[index3] = sprintf(
      "%s card or history results of %.0f percent modified for recall bias to %.0f percent based on 1st dose card or history coverage of %.0f percent, 1st dose card only coverage of %.0f percent and 3rd dose card only coverage of %.0f percent. ", 
      Svy.Title[index3], tround(Svy.Ana[index3]), H3Adj[index3], tround(Svy.CH1[index1]), tround(Svy.C1[index1]), tround(Svy.C3[index3]))
    Svy.Ana[index3] = H3Adj[index3]
  }
  
  # Some surveys are ignored by the working group
  ignore = Decisions[Decisions$Dec == "ignoreSurvey" & !is.na(Decisions$Id), ]
  Svy.Acc = Svy.Ana
  Svy.Acc[cbind(ignore$Y, ignore$V, ignore$Id)] = NA
  # commented out, keep information
  # Svy.Info[cbind(ignore$Y, ignore$V, ignore$Id)] = NA
  
  # Some surveys are ignored by the working group (by year and vaccine, no Id)
  ignore = Decisions[Decisions$Dec == "ignoreSurvey" & is.na(Decisions$Id), ]
  if(nrow(ignore))
    for(i in 1:nrow(ignore))
    {
      Svy.Acc[ignore$Y[i], ignore$V[i], ] = NA
      # commented out, keep information
      # Svy.Info[ignore$Y[i], ignore$V[i], ] = NA
    }
  
  # Workgroup decision: keep
  accept   = Decisions[Decisions$Dec == "acceptSurvey" & !is.na(Decisions$Id), ]
  
  cnf      = Survey$Info.confirm == "card or history"
  age      = Survey$Info.age %in% c("12-23 m", "18-29 m", "15-26 m", "24-35 m")
  size     = Survey$Info.ss < 300
  
  if(any(cnf & age & size))
  {
    ignore   = Survey[cnf & age & size, ]
    
    # remove "ignore" if it is already in accept
    if(nrow(accept))
    {
      acc_ign = rbind(accept[, c("Y", "V", "Id")], ignore[, c("Y", "V", "Id")]) 
      ignore = unique(acc_ign)[-(1:nrow(accept)), , drop=FALSE]
    }
    Svy.Acc[cbind(ignore$Y, ignore$V, ignore$Id)] = NA
  }

  # % Survey information for given year. Multiple surveys are averaged.
  # survey(C, V, Y, Expl, Coverage) :-
  #     findall(Cov, survey_accepted(C, V, Y, _, Cov), [H | T]),
  #     length([H | T], N),
  #     sum_list([H | T], Sum),
  #     Coverage is round(Sum / N),
  #     concat_atom(['Survey evidence of ', Coverage, ' percent based on ',
  #         N, ' survey(s). '], Expl).
  
  Svy.Cov = tround(apply(Svy.Acc, c(1, 2), mean, na.rm=TRUE))
  Svy.Expl[] = sprintf(
    "Survey evidence of %g percent based on %i survey(s). ", 
    Svy.Cov, apply(!is.na(Svy.Acc), c(1, 2), sum))
}


# 6. Determine coverage value at anchor points defined as years with multiple
#    data points (reported | survey | wgd).

Anchor.Rule = YV_char
Anchor.Info = YV_char
Anchor.Cov = YV_int

# % Survey results challenge reported
# anchor(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, _, _Cov0),
#     survey(C, V, Y, Expl0, Survey),
#     % survey_reported_threshold(Threshold),
#     % abs(Cov0 - Survey) > Threshold,
#     !,
#     Rule = 'S: AP',
#     concat_atom(['Survey evidence does not support reported data. Estimate based on survey results. ',
#     Expl0, ' '], Expl),
#     Coverage = Survey.

index = which(abs(Svy.Cov - TS.Cov) > svy.thrs, arr.ind=TRUE)
Anchor.Rule[index] = "S: AP"
Anchor.Info[index] = sprintf(
  "Survey evidence does not support reported data. Estimate based on survey results. %s ",
  Svy.Expl[index])
Anchor.Cov[index] = Svy.Cov[index]

# % Survey results support reported
# anchor(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, Source, Cov0),
#     survey(C, V, Y, Expl0, Survey),
#     survey_reported_threshold(Threshold),
#     abs(Cov0 - Survey) =< Threshold,
#     !,
#     Rule = 'R: AP',
#     member(Source-Expl1,
#       [ gov-'Estimate informed by reported data supported by survey. ',
#         admin-'Estimate informed by reported administrative data supported by survey. ',
#         interpolated-'Estimate informed by interpolation between reported data supported by survey. ',
#         extrapolated-'Estimate based on extrapolation from data reported by national government supported by survey. '
#       ]),
#     concat_atom([Expl1, Expl0], Expl),
#     Coverage = Cov0.

info = c(
  gov="Estimate informed by reported data supported by survey. ",
  admin="Estimate informed by reported administrative data supported by survey. ",
  interpolated="Estimate informed by interpolation between reported data supported by survey. ",
  extrapolated="Estimate based on extrapolation from data reported by national government supported by survey. ")

index = which(abs(Svy.Cov - TS.Cov) <= svy.thrs, arr.ind=TRUE)
Anchor.Rule[index] = "R: AP"
Anchor.Cov[index] = TS.Cov[index]
Anchor.Info[index] = sprintf("%s%s", info[TS.Src[index]], Svy.Expl[index])

# % Reported value "anchored" by working group
# anchor(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, _, Cov0),
#     decision(C, V, Y, assignAnchor, Expl0, _, Cov0), % same Cov0
#     !,
#     Rule = 'R: AP',
#     Expl = Expl0,
#     Coverage = Cov0.

index = Decisions[Decisions$Dec == "assignAnchor", ]
equal = which(TS.Cov[cbind(index$Y, index$V)] == index$Cov, arr.ind=TRUE)
index = index[equal, ]
Anchor.Rule[cbind(index$Y, index$V)] = "R: AP"
Anchor.Info[cbind(index$Y, index$V)] = index$Info
Anchor.Cov[cbind(index$Y, index$V)] = index$Cov

# % Working group assigns anchor point value.
# anchor(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, _, _Cov0),
#     decision(C, V, Y, assignAnchor, Expl0, _, Assigned),
#     % Cov0 \= Assigned,
#     !,
#     Rule = 'W: AP',
#     concat_atom(['Estimate of ', Assigned, ' percent assigned by working group. ',
#         Expl0], Expl),
#     Coverage = Assigned.

index = Decisions[Decisions$Dec == "assignAnchor", ]
neq = which(TS.Cov[cbind(index$Y, index$V)] != index$Cov, arr.ind=TRUE)
index = index[neq, ]
Anchor.Rule[cbind(index$Y, index$V)] = "W: AP"
Anchor.Cov[cbind(index$Y, index$V)] = index$Cov
Anchor.Info[cbind(index$Y, index$V)] = sprintf(
  "Estimate of %g percent assigned by working group. %s", index$Cov, index$Info)

# % 7. Estimate coverage by distinguishing different cases
# %
# % * Estimate at anchor point
# % * Estimate between anchor points
# % * Estimate before anchor point
# % * Estimate after anchor point
# % * No anchor points
#
# % No anchor points for any year, use reported
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, Source, Cov0),
#     !,
#     Rule = 'R:',
#     member(Source-Expl,
#       [ gov-'Estimate informed by reported data. ',
#         admin-'Estimate informed by reported administrative data. ',
#         interpolated-'Estimate informed by interpolation between reported data. ',
#         extrapolated-'Estimate informed by extrapolation from reported data. '
#       ]),
#     Coverage = Cov0.

info = c(
  gov="Estimate informed by reported data. ",
  admin="Estimate informed by reported administrative data. ",
  interpolated="Estimate informed by interpolation between reported data. ",
  extrapolated="Estimate informed by extrapolation from reported data. ")
Cov = TS.Cov

Info = YV_char
Info[] = info[TS.Src]

Rule = YV_char
Rule[!is.na(TS.Cov)] = "R:"

# % Before earliest/after latest anchor (not of type reported): calibrated
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, _, Cov0),
#     (   succ_anchor(C, V, Y, Anchor, AnchorRule, AnchorCov),
#         AnchorRule \= 'R: AP'
#     ;   prec_anchor(C, V, Y, Anchor, AnchorRule, AnchorCov),
#         AnchorRule \= 'R: AP'
#     ), !,
#     Rule = 'C:',
#     concat_atom(['Reported data calibrated to ', Anchor, ' levels. '], Expl),
#     reported_time_series(C, V, Anchor, _, ReportedAtAnchor),
#     Adj is AnchorCov - ReportedAtAnchor,
#     Coverage is round(Cov0 + Adj).

# Search for preceding anchor
index = !is.na(TS.Cov) & is.na(Anchor.Cov)
Prec.Rule = apply(Anchor.Rule, 2, FUN=na.locf, na.rm=FALSE)
index = index & !is.na(Prec.Rule) & Prec.Rule != "R: AP"

Prec.Cov = apply(Anchor.Cov, 2, FUN=na.locf, na.rm=FALSE)

Prec.Year = YV_char
Prec.Year[] = Yn
Prec.Year[index] = NA
Prec.Year = apply(Prec.Year, 2, FUN=na.locf, na.rm=FALSE)

# Rule = YV_char
Rule[index] = "C:"

# Info = YV_char
Info[index] = sprintf("Reported data calibrated to %s levels. ", 
    Prec.Year[index])

yv = expand.grid(Y=Yn, V=Vn, stringsAsFactors=FALSE)
Adj = YV_int
Adj[] = Anchor.Cov[cbind(c(Prec.Year), yv$V)] - TS.Cov[cbind(c(Prec.Year), yv$V)]
# Cov = YV_int
Cov[index] = TS.Cov[index] + Adj[index]

# Search for next anchor
index = !is.na(TS.Cov) & is.na(Anchor.Cov)
Succ.Rule = apply(Anchor.Rule, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)
index = index & !is.na(Succ.Rule) & Succ.Rule != "R: AP"

Succ.Cov = apply(Anchor.Cov, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)

Succ.Year = YV_char
Succ.Year[] = Yn
Succ.Year[index] = NA
Succ.Year = apply(Succ.Year, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)

Rule[index] = "C:"

Info[index] = sprintf("Reported data calibrated to %s levels. ", 
    Succ.Year[index])

yv = expand.grid(Y=Yn, V=Vn, stringsAsFactors=FALSE)
Adj = Anchor.Cov[cbind(c(Succ.Year), yv$V)] - TS.Cov[cbind(c(Succ.Year), yv$V)]
Cov[index] = TS.Cov[index] + Adj[index]

# % Before earliest/after latest anchor (of type reported)
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, Source, Cov0),
#     (   succ_anchor(C, V, Y, _Anchor, AnchorRule, _AnchorCov),
#         AnchorRule = 'R: AP'
#     ;   prec_anchor(C, V, Y, _Anchor, AnchorRule, _AnchorCov),
#         AnchorRule = 'R: AP'
#     ), !,
#     Rule = 'R:',
#     member(Source-Expl,
#       [ gov-'Estimate informed by reported data. ',
#         admin-'Estimate informed by reported administrative data. ',
#         interpolated-'Estimate informed by interpolation between reported data. ',
#         extrapolated-'Estimate based on extrapolation from data reported by national government. '
#       ]),
#     Coverage = Cov0.

# Search for preceding anchor
index = !is.na(TS.Cov) & is.na(Anchor.Cov)
Prec.Rule = apply(Anchor.Rule, 2, FUN=na.locf, na.rm=FALSE)
index = index & !is.na(Prec.Rule) & Prec.Rule == "R: AP"

info = c(gov="Estimate informed by reported data. ",
    admin="Estimate informed by reported administrative data. ",
    interpolated="Estimate informed by interpolation between reported data. ",
    extrapolated="Estimate based on extrapolation from data reported by national government. ")

Rule[index] = "R:"
Info[index] = info[TS.Src[index]]
Cov[index] = TS.Cov[index]

# Search for next anchor
index = !is.na(TS.Cov) & is.na(Anchor.Cov)
Succ.Rule = apply(Anchor.Rule, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)
index = index & !is.na(Succ.Rule) & Succ.Rule == "R: AP"

Rule[index] = "R:"
Info[index] = info[TS.Src[index]]
Cov[index] = TS.Cov[index]

# % Between other anchor points (not both of type "reported"): calibrate
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, _Source, _Cov0),
#     prec_anchor(C, V, Y, Prec, PrecRule, PrecCov),
#     succ_anchor(C, V, Y, Succ, SuccRule, SuccCov),
#     ( PrecRule \= 'R: AP' ; SuccRule \= 'R: AP' ),
#     !,
#     Rule = 'C:',
#     concat_atom(['Reported data calibrated to ', Prec,
#         ' and ', Succ, ' levels. '], Expl),
#     reported_time_series(C, V, Prec, _, PrecRep),
#     reported_time_series(C, V, Succ, _, SuccRep),
#     interpolate(Prec, PrecRep, Succ, SuccRep, Y, RepInterp),
#     interpolate(Prec, PrecCov, Succ, SuccCov, Y, AnchInterp),
#     Adj is AnchInterp - RepInterp,
#     Coverage is round(Reported + Adj).

index = !is.na(TS.Cov) & is.na(Anchor.Cov)
Prec.Rule = apply(Anchor.Rule, 2, FUN=na.locf, na.rm=FALSE)
Succ.Rule = apply(Anchor.Rule, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)
index = index & !is.na(Prec.Rule) & !is.na(Succ.Rule)
index = index & (Prec.Rule != "R: AP" | Succ.Rule != "R: AP")

Prec.Cov = apply(Anchor.Cov, 2, FUN=na.locf, na.rm=FALSE)
Succ.Cov = apply(Anchor.Cov, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)

Rule[index] = "C:"

Prec.Year = YV_char
Prec.Year[] = Yn
Prec.Year[index] = NA
Prec.Year = apply(Prec.Year, 2, FUN=na.locf, na.rm=FALSE)

Succ.Year = YV_char
Succ.Year[] = Yn
Succ.Year[index] = NA
Succ.Year = apply(Succ.Year, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)

Info[index] = sprintf("Reported data calibrated to %s and %s levels. ", 
    Prec.Year[index], Succ.Year[index])

# interpolate Anchor.Cov for year without anchor
Itp1.Cov = apply(Anchor.Cov, 2, FUN=na.approx, na.rm=FALSE)

# MG, discuss: this shouldn't be rounded
Itp1.Cov[] = tround(Itp1.Cov)
rownames(Itp1.Cov) = Yn

# interpolate TS.Cov for year without anchor
Itp2.Cov = TS.Cov
Itp2.Cov[index] = NA
Itp2.Cov = apply(Itp2.Cov, 2, FUN=na.approx, na.rm=FALSE)
rownames(Itp2.Cov) = Yn
Itp2.Cov[] = tround(Itp2.Cov)

# yv = expand.grid(Y=Yn, V=Vn, stringsAsFactors=FALSE)
Adj = Itp1.Cov - Itp2.Cov
Cov[index] = tround(TS.Cov[index] + Adj[index])

# % Between anchor points: between two reported anchors
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     reported_time_series(C, V, Y, Source, Cov0),
#     prec_anchor(C, V, Y, _Prec, PrecRule, _),
#     PrecRule = 'R: AP',
#     succ_anchor(C, V, Y, _Succ, SuccRule, _),
#     SuccRule = 'R: AP',
#     !,
#     Rule = 'R:',
#     member(Source-Expl,
#       [ gov-'Estimate informed by reported data. ',
#         admin-'Estimate informed by reported administrative data. ',
#         interpolated-'Estimate informed by interpolation between reported data. '
#       ]),
#     Coverage = Cov0.

index = !is.na(TS.Cov) & is.na(Anchor.Cov)
Prec.Rule = apply(Anchor.Rule, 2, FUN=na.locf, na.rm=FALSE)
Succ.Rule = apply(Anchor.Rule, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)
index = index & !is.na(Prec.Rule) & !is.na(Succ.Rule)
index = index & Prec.Rule == "R: AP" & Succ.Rule == "R: AP"

index1 = index & TS.Src %in% c("gov", "admin", "interpolated")
info = c(gov="Estimate informed by reported data. ",
    admin="Estimate informed by reported administrative data. ",
    interpolated="Estimate informed by interpolation between reported data. ")

Rule[index1] = "R:"
Info[index1] = info[TS.Src[index1]]
Cov[index1] = TS.Cov[index1]

# MG, discuss: if it is extrapolated, the rule fails and NA is returned. See the
# corresponding Prolog rule.
index2 = index & TS.Src == "extrapolated"
Rule[index2] = NA
Info[index2] = NA
Cov[index2] = NA

# % Between anchor points: interpolation forced by working group
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     decision(C, V, Y, interpolate, Expl0, _, _),
#     prec_anchor(C, V, Y, Prec, _, PrecCov),
#     succ_anchor(C, V, Y, Succ, _, SuccCov),
#     !,
#     Rule = 'W-I:',
#     concat_atom(['Estimate informed by interpolation between ', Prec,
#         ' and ', Succ, ' levels. ', Expl0], Expl),
#     interpolate(Prec, PrecCov, Succ, SuccCov, Y, Coverage).

index = Decisions[Decisions$Dec == "interpolate", ]

Prec.Rule = apply(Anchor.Rule, 2, FUN=na.locf, na.rm=FALSE)
Prec.Cov = apply(Anchor.Cov, 2, FUN=na.locf, na.rm=FALSE)
Succ.Rule = apply(Anchor.Rule, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)
Succ.Cov = apply(Anchor.Cov, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)

Prec.Year = YV_char
Prec.Year[] = Yn
Prec.Year[is.na(Anchor.Cov)] = NA
# Prec.Year[cbind(index$Y, index$V)] = NA
Prec.Year = apply(Prec.Year, 2, FUN=na.locf, na.rm=FALSE)

Succ.Year = YV_char
Succ.Year[] = Yn
Succ.Year[is.na(Anchor.Cov)] = NA
# Succ.Year[cbind(index$Y, index$V)] = NA
Succ.Year = apply(Succ.Year, 2, FUN=na.locf, fromLast=TRUE, na.rm=FALSE)

WI = YV_char
WI[cbind(index$Y, index$V)] = index$Info
WI[WI == "NA"] = ""
WI[is.na(Prec.Rule) | is.na(Succ.Rule)] = NA

index = !is.na(WI)
Rule[index] = "W-I:"
Info[index] = sprintf(
  "Estimate informed by interpolation between %s and %s levels. %s",
  Prec.Year[index], Succ.Year[index], WI[index])

Itp1.Cov = apply(Anchor.Cov, 2, FUN=na.approx, na.rm=FALSE)
rownames(Itp1.Cov) = rownames(Anchor.Cov)
Cov[index] = tround(Itp1.Cov[index])

# % At anchor points
# wuenic_II(C, V, Y, Rule, Expl, Coverage) :-
#     anchor(C, V, Y, Rule0, Expl0, Cov0),
#     !,
#     Rule = Rule0,
#     Expl = Expl0,
#     Coverage = Cov0.

index = !is.na(Anchor.Cov)
Rule[index] = Anchor.Rule[index]
Info[index] = Anchor.Info[index]
Cov[index] = Anchor.Cov[index]

# % Obtain estimates for vaccine coverage. At Level 1, check for working group
# % decisions and work around obvious inconsistencies.
#
# % Estimate for RCV1 where RCV1 given at MCV1
# wuenic_I(C, rcv1, Y, Rule, Expl, Coverage) :-
#   estimate_required(C, rcv1, Y, _, _),
#   !,
#   wuenic_II(C, mcv1, Y, Rule, _, Coverage),
#   Expl = 'Estimate based on estimated MCV1. '.

index = Ereq[, "rcv1"]
Rule[index, "rcv1"] = Rule[index, "mcv1"]
Info[index, "rcv1"] = "Estimate based on estimated MCV1. "
Cov[index, "rcv1"] = Cov[index, "mcv1"]

# % Estimate for RCV1 where RCV1 given at MCV2
# wuenic_I(C, rcv1, Y, Rule, Expl, Coverage) :-
#     estimate_required(C, rcv1, Y, _, FirstRubellaDose),
#     firstRubellaAtSecondMCV(C, rcv1, Y, FirstRubellaDose),
#     !,
#     wuenic_II(C, mcv2, Y, Rule, _, Coverage),
#     Expl = 'First dose of rubella vaccine given with second dose of measles containing vaccine. Estimate based on MCV2 estimate'.

index = Rubella[, "rcv1"] == firstRubellaAtSecondMCV[, "rcv1"]
index = as.character(Yn[which(index)])
if(length(index))
{
  Rule[cbind(index, "rcv1")] = Rule[cbind(index, Rubella[, "rcv1"][index])]
  Info[cbind(index, "rcv1")] = "First dose of rubella vaccine given with second dose of measles containing vaccine. Estimate based on MCV2 estimate"
  Cov[cbind(index, "rcv1")] = Cov[cbind(index, Rubella[, "rcv1"][index])]
}

# % If DTP1 not reported: estimate using equation
# % If DTP3 > DTP1 (which is impossible), estimate coverage using equation
# wuenic_I(C, dtp1, Y, Rule, Expl, Coverage) :-
#     wuenic_II(C, dtp3, Y, _, _, DTP3),
#     !,
#     Rule = 'RMF:',
#     concat_atom(['Estimate based on DTP3 coverage of ', DTP3, '. '], Expl),
#     Coverage is round(-0.0058 * DTP3 * DTP3 + 1.3912 * DTP3 + 18.258).

index = !is.na(Cov[, "dtp1"]) & !is.na(Cov[, "dtp3"]) & Cov[, "dtp3"] > Cov[, "dtp1"]
index = index | (is.na(Cov[, "dtp1"]) & !is.na(Cov[, "dtp3"]))
Rule[index, "dtp1"] = "RMF:"
Info[index, "dtp1"] = sprintf("Estimate based on DTP3 coverage of %i. ",
  Cov[index, "dtp3"])
Cov[index, "dtp1"] =
  round(-0.0058 * Cov[index, "dtp3"]^2 + 1.3912 * Cov[index, "dtp3"] + 18.258)

# % Assigned by working group
# wuenic_I(C, V, Y, Rule, Expl, Coverage) :-
#     decision(C, V, Y, assignWUENIC, Expl0, _, Cov0),
#     !,
#     Rule = 'W:',
#     Expl = Expl0,
#     Coverage = Cov0.

index = Decisions[Decisions$Dec == "assignWUENIC", ]
Rule[cbind(index$Y, index$V)] = "W:"
Info[cbind(index$Y, index$V)] = index$Info
Cov[cbind(index$Y, index$V)] = index$Cov

# 9. Force between 0 and 99%
#
# bound_0_100(X, Y) :-
#   Y is max(0, min(99, round(X))).
Bounded = Cov
Bounded[] = pmax(0, pmin(99, round(Bounded)))

# % Confidence in reported coverage if it is an anchor point. Otherwise,
# % no confidence
# conf_reported(C, V, Y, Support) :-
#   reported(C, V, Y, _, _),
#   wuenic_I(C, V, Y, Rule, _, _),
#   !,
#   (   member(Rule, ['R:', 'R: AP'])
#   ->  Support = 'R+'
#   ;   Support = 'R-'
#   ).

Conf.Rep = YV_char
index = !is.na(Rep.Cov) & Rule %in% c("R:", "R: AP")
Conf.Rep[index] = "R+"
index = !is.na(Rep.Cov) & !(Rule %in% c("R:", "R: AP"))
Conf.Rep[index] = "R-"

# % No confidence in surveys if _any_ survey deviates too much from WUENIC
# % coverage
# conf_survey(C, V, Y, Support) :-
#   wuenic_I(C, V, Y, _, _, Cov0),
#   estimate_required(C, V, Year, _, _),
#   survey(C, V, Year, _, Coverage),
#   confidence_survey_scope(Scope),
#   abs(Y - Year) =< Scope,
#   confidence_survey_threshold(Threshold),
#   abs(Cov0 - Coverage) > Threshold,
#   !,
#   Support = 'S-'.
#
# % Confidence only if all surveys are consistent with WUENIC coverage
# conf_survey(C, V, Y, Support) :-
#   wuenic_I(C, V, Y, _, _, _Cov0),
#   estimate_required(C, V, Year, _, _),
#   survey(C, V, Year, _, _Coverage),
#   confidence_survey_scope(Scope),
#   abs(Y - Year) =< Scope,
#   !,
#   Support = 'S+'.

# MG, discuss: It seems as if surveys are ignored if there is no estimate
# required. Unclear if intended, since the scope of surveys spans +/- 2 years.
# Example: arm/pcv3 2015 (Svy.Min from 2013)
index = which(!Ereq, arr.ind=TRUE)

Svy.Min = Svy.Cov
Svy.Min[cbind(index[, "Y"], index[, "V"])] = NA # here
Svy.Min = suppressWarnings({rollapply(Svy.Min, width=1 + 2*svy.scope, partial=1, FUN=min, na.rm=TRUE)})
rownames(Svy.Min) = rownames(Svy.Cov)

Svy.Max = Svy.Cov
Svy.Max[cbind(index[, "Y"], index[, "V"])] = NA # here
Svy.Max = suppressWarnings({rollapply(Svy.Max, width=1 + 2*svy.scope, partial=1, FUN=max, na.rm=TRUE)})
rownames(Svy.Max) = rownames(Svy.Cov)

Conf.Svy = YV_char
index = which(!is.na(Cov) & Ereq & is.finite(Svy.Min) & is.finite(Svy.Max))
Conf.Svy[index] = "S+"

index = which(!is.na(Cov) & Ereq & Cov - Svy.Min > svy.thrs)
Conf.Svy[index] = "S-"

index = which(!is.na(Cov) & Ereq & Svy.Max - Cov > svy.thrs)
Conf.Svy[index] = "S-"

# % Births used for bcg and hepb birth dose
# denominator(C, V, Y, Coverage) :-
#   member(V, [bcg, hepbb]),
#   !,
#   vaccinated0(C, V, Y, Vaccinated),
#   births_UNPD(C, Y, Births),
#   Coverage is Vaccinated / Births * 100.
#
# % Surviving infants for remaining vaccines
# denominator(C, V, Y, Coverage) :-
#   vaccinated0(C, V, Y, Vaccinated),
#   si_UNPD(C, Y, SI),
#   Coverage is Vaccinated / SI * 100.

Den = YV_real
index = c("bcg", "hepbb")
Den[, index] = Vaccinated[, index] / Births * 100

index = setdiff(Vn, c("bcg", "hepbb"))
Den[, index] = Vaccinated[, index] / Surviving * 100

# % Recalculate coverage using reported number of children vaccinated and
# % births and surviving infants from UNPD estimates.
# conf_denominator(C, V, Y, Support) :-
#   vaccinated0(C, V, Y, _),
#   births_UNPD(C, Y, _),
#   si_UNPD(C, Y, _),
#   wuenic_I(C, V, Y, _Rule, _Expl, Cov0),
#   denominator(C, V, Y, Coverage),
#   !,
#   confidence_UNPD_threshold(Threshold),
#   (   abs(Coverage - Cov0) < Threshold % MG: < inconsistent with conf_survey
#   ->  Support = 'D+'
#   ;   Support = 'D-'
#   ).

Conf.Den = YV_char
index = which(abs(Cov - Den) < unpd.thrs)
Conf.Den[index] = "D+"
index = which(abs(Cov - Den) >= unpd.thrs)
Conf.Den[index] = "D-"

# 10. Confidence depends on converging evidence from the different sources.
# % 1 star = low confidence, ..., 3 stars = high confidence
# %
# % Low confidence
# confidence(_C, _V, _Y, Expl, Grade) :-
#   !,
#   Expl = 'GoC=No accepted empirical data',
#   Grade = 1.

GoC = YV_int
GoC.Expl = YV_char
GoC[] = 1
GoC.Expl[] = "GoC=No accepted empirical data"

# % Confidence in one or two sources, two stars
# confidence(C, V, Y, Expl, Grade) :-
#   conf_reported(C, V, Y, 'R+'),
#   conf_survey(C, V, Y, 'S+'),
#   !,
#   Expl = 'GoC=R+ S+',
#   Grade = 2.
#
# confidence(C, V, Y, Expl, Grade) :-
#   conf_survey(C, V, Y, 'S+'),
#   conf_denominator(C, V, Y, 'D+'),
#   !,
#   Expl = 'GoC=S+ D+',
#   Grade = 2.
#
# confidence(C, V, Y, Expl, Grade) :-
#   conf_reported(C, V, Y, 'R+'),
#   conf_denominator(C, V, Y, 'D+'),
#   !,
#   Expl = 'GoC=R+ D+',
#   Grade = 2.
#
# confidence(C, V, Y, Expl, Grade) :-
#   conf_reported(C, V, Y, 'R+'),
#   !,
#   Expl = 'GoC=R+',
#   Grade = 2.
#
# confidence(C, V, Y, Expl, Grade) :-
#   conf_survey(C, V, Y, 'S+'),
#   !,
#   Expl = 'GoC=S+',
#   Grade = 2.
#
# confidence(C, V, Y, Expl, Grade) :-
#   conf_denominator(C, V, Y, 'D+'),
#   !,
#   Expl = 'GoC=D+',
#   Grade = 2.

index = which(Conf.Rep == "R+")
GoC[index] = 2
GoC.Expl[index] = "GoC=R+"

index = which(Conf.Svy == "S+")
GoC[index] = 2
GoC.Expl[index] = "GoC=S+"

index = which(Conf.Den == "D+")
GoC[index] = 2
GoC.Expl[index] = "GoC=D+"

index = which(Conf.Rep == "R+" & Conf.Svy == "S+")
GoC[index] = 2
GoC.Expl[index] = "GoC=R+ S+"

index = which(Conf.Svy == "S+" & Conf.Den == "D+")
GoC[index] = 2
GoC.Expl[index] = "GoC=S+ D+"

index = which(Conf.Rep == "R+" & Conf.Den == "D+")
GoC[index] = 2
GoC.Expl[index] = "GoC=R+ D+"

# % Check if any source is challenged
# challenge(C, V, Y, 'R-') :-
#   conf_reported(C, V, Y, 'R-').
#
# challenge(C, V, Y, 'S-') :-
#   conf_survey(C, V, Y, 'S-').
#
# challenge(C, V, Y, 'D-') :-
#   conf_denominator(C, V, Y, 'D-').
#
# % Todo list from V3
# %
# % 1. Simplify previous rule to a check for an anchor point
# %
# % supporting_survey_in_scope(C, V, Y, Rule) :-
# %     survey(C, V, Y, _, _),
# %     wuenic_I(C, V, Y, 'S: AP', _, _).
# %
# % 2. Rewrite rule to look at relationship between estimate rule and
# %    surveys in scope rule. For example, take randomness into account,
# %    e.g. probability for inconsistent results increases with the number of
# %    surveys and decreases with the sample size of the surveys.

Chall = YV_char
Chall[] = ""

index = which(Conf.Den == "D-")
Chall[index] = "D-"

index = which(Conf.Rep == "R-")
Chall[index] = sprintf("%sR-", Chall[index])

index = which(Conf.Svy == "S-")
Chall[index] = sprintf("%sS-", Chall[index])

# % If any estimate has been challenged, confidence is low
# confidence(C, V, Y, Expl, Grade) :-
#   setof(Expl0, challenge(C, V, Y, Expl0), List),
#   !,
#   concat_atom(['Estimate challenged by: ' | List], Expl),
#   Grade = 1.

index = Chall != ""
GoC[index] = 1
GoC.Expl[index] = sprintf("Estimate challenged by: %s", Chall[index])

# % Confidence in both reported, surveys, and sold vaccines
# confidence(C, V, Y, Expl, Grade) :-
#   conf_reported(C, V, Y, 'R+'),
#   conf_survey(C, V, Y, 'S+'),
#   conf_denominator(C, V, Y, 'D+'),
#   !,
#   Expl = 'GoC=R+ S+ D+',
#   Grade = 3.

index = which(Conf.Rep == "R+" & Conf.Svy == "S+" & Conf.Den =="D+")
GoC[index] = 3
GoC.Expl[index] = "GoC=R+ S+ D+"

# % Confidence rated by working group
# confidence(C, V, Y, Expl, Grade) :-
#   decision(C, V, Y, assignGoC, Expl0, _, Grade0),
#   !,
#   concat_atom(['GoC=Assigned by working group. ', Expl0], Expl),
#   Grade = Grade0.

index = Decisions[Decisions$Dec == "assignGoC", ]
GoC[cbind(index$Y, index$V)] = index$Cov # Cov is GoC here
GoC.Expl[cbind(index$Y, index$V)] = sprintf("GoC=Assigned by working group. %s",
    index$Info)

# % Copy rcv1 from mcv2
# %
# % MG, discuss: this differs from the rule in wuenic_I
# confidence(C, rcv1, Y, Expl, Grade) :-
#   estimate_required(C, rcv1, Y, _, mcv2),
#   !,
#   confidence(C, mcv2, Y, Expl, Grade).

index = which(Ereq[, "rcv1"] & Rubella[, "rcv1"] == "mcv2")
GoC[index, "rcv1"] = GoC[index, "mcv2"]
GoC.Expl[index, "rcv1"] = GoC.Expl[index, "mcv2"]

# % Copy rcv1 from mcv1
# confidence(C, rcv1, Y, Expl, Grade) :-
#   estimate_required(C, rcv1, Y, _, na),
#   !,
#   confidence(C, mcv1, Y, Expl, Grade).

index = Ereq[, "rcv1"] & is.na(Rubella[, "rcv1"])
GoC[index, "rcv1"] = GoC[index, "mcv1"]
GoC.Expl[index, "rcv1"] = GoC.Expl[index, "mcv1"]

# % Flag modifications in the program code that change the coverage
# % estimates
# change_from_previous(C, V, Y, Coverage, Change) :-
#     legacy(C, V, Y, Legacy),
#     Legacy \= Coverage,
#     !,
#     concat_atom(['Estimate of ', Coverage,
#         ' percent changed from previous revision value of ',
#         Legacy,' percent. '], Change).
#
# change_from_previous(_C, _V, _Y, _, '').

Change = YV_char
Change[] = ""
index = which(Bounded != Legacy)
Change[index] = sprintf(
  "Estimate of %i percent changed from previous revision value of %i percent. ",
  Bounded[index], Legacy[index])

# % Collect explanations in natural language terms
# collect_explanations(C, V, Y, Explanations) :-
#     findall(Expl, explanation(C, V, Y, Expl), Explanations).

Expl = YV_char
Expl[] = ""

# explanation(C, V, Y, Expl) :-
#     survey_reason_to_exclude(C, V, Y, _, Expl).
#
# % Reasons to exclude a survey include:
# %    Sample size < 300,
# %    The working group decides to exclude the survey.
# %
# % used for explanation/4
# survey_reason_to_exclude(C, V, Y, ID, Expl) :-
#     survey_for_analysis(C, V, Y, ID, Description, _),
#     not(decision(C, V, Y, acceptSurvey, Expl, ID, _)),
#     member(ss:Size, Description),
#     Size < 300,
#     concat_atom(['Survey results ignored. Sample size ', Size,
#         ' less than 300. '], Expl).

cnf    = Survey$Info.confirm == "card or history"
age    = Survey$Info.age %in% c("12-23 m", "18-29 m", "15-26 m", "24-35 m")
# accept = Survey$Id %in% Decisions$Id[Decisions$Dec == "acceptSurvey"]
size   = Survey$Info.ss < 300

if(any(cnf & age & size))
{
  index  = Survey[cnf & age & size, ]
  accept   = Decisions[Decisions$Dec == "acceptSurvey" & !is.na(Decisions$Id), ]
  
  # MG, discuss: bgd, bcg/1999: twice the same message, without Survey ID
  for(i in 1:nrow(index))
  {
    if(!any(index$V[i] == accept$V & index$Y[i] == accept$Y & index$Id[i] == accept$Id))
      Expl[index$Yn[i], index$V[i]] = sprintf(
        "%sSurvey results ignored. Sample size %i less than 300. ", 
        Expl[index$Yn[i], index$V[i]], index$Info.ss[i])
  }
}

# % V4: Keeps explanations for the same survey together
# survey_reason_to_exclude(C, V, Y, ID, Expl) :-
#     survey_for_analysis(C, V, Y, ID, Description, _),
#     (   decision(C, V, Y, ignoreSurvey, Expl0, ID, _)
#     ;   decision(C, V, Y, ignoreSurvey, Expl0, na, _)
#     ),
#     member(title:Title, Description),
#     concat_atom([Title, ' results ignored by working group. ', Expl0],
#     Expl).

# Avoid duplicate "ignored by working group"
Ignored = array("", dim=c(length(Yn), length(Vn), length(Dn)), 
                dimnames=list(Y=Yn, V=Vn, Id=Dn))

# Some surveys are ignored by the working group (by year and vaccine, no Id)
ignore = Decisions[Decisions$Dec == "ignoreSurvey" & is.na(Decisions$Id), ]
if(nrow(ignore))
  for(i in 1:nrow(ignore))
  {
    Ids = Dn[!is.na(Svy.Ana[ignore$Y[i], ignore$V[i], , drop=FALSE])]
    if(length(Ids))
      for(j in Ids)
        Ignored[ignore$Y[i], ignore$V[i], j] = sprintf(
          "%s%s results ignored by working group.%s ",
          Ignored[ignore$Y[i], ignore$V[i], j],
          Survey$Info.title[Survey$Id == j][1], ignore$Info[i])
  }

# Some surveys are ignored by the working group
ignore = Decisions[Decisions$Dec == "ignoreSurvey" & !is.na(Decisions$Id), ]
if(nrow(ignore))
  for(i in 1:nrow(ignore))
    if(!is.na(Svy.Ana[ignore$Y[i], ignore$V[i], ignore$Id[i]]))
      Ignored[ignore$Y[i], ignore$V[i], ignore$Id[i]] = sprintf(
        "%s%s results ignored by working group.%s ",
        Ignored[ignore$Y[i], ignore$V[i], ignore$Id[i]],
        Survey$Info.title[Survey$Id == ignore$Id[i]][1], ignore$Info[i])

Ignored = apply(Ignored, c(1, 2), paste, collapse="")
Expl[] = sprintf("%s%s", Expl, Ignored)

# explanation(C, V, Y, Expl) :-
#     survey_modified(C, V, Y, _, Expl, _).

index = which(!is.na(Svy.Info), arr.ind=TRUE)
if(nrow(index))
  for(i in 1:nrow(index))
  {
    Expl[index[i, "Y"], index[i, "V"]] = sprintf("%s%s",
      Expl[index[i, "Y"], index[i, "V"]],
      Svy.Info[index[i, "Y"], index[i, "V"], index[i, "Id"]])
  }

# explanation(C, V, Y, Expl) :-
#     reported_reason_to_exclude(C, V, Y, _, Expl).
#
# % Reasons to exclude reported data are:
# %  1. Working group decision.
# %  2. Coverage > 100%
# %  3. Inconsistent temporal changes (sawtooth or sudden change most
# %     recent year)
# reported_reason_to_exclude(C, V, Y, Reason, Expl) :-
#     reported(C, V, Y, _, _),
#     decision(C, V, Y, ignoreReported, Expl0, _, _),
#     Reason = wdg,
#     concat_atom(['Reported data excluded. ', Expl0], Expl).

ignore = Decisions[Decisions$Dec == "ignoreReported", ]
if(nrow(ignore))
  for(i in 1:nrow(ignore))
  {
    if(!is.na(Rep.Cov[ignore$Y[i], ignore$V[i]]))
      Expl[ignore$Y[i], ignore$V[i]] = sprintf("%sReported data excluded. %s",
        Expl[ignore$Y[i], ignore$V[i]], ignore$Info[i])
  }


# reported_reason_to_exclude(C, V, Y, Reason, Expl) :-
#     reported(C, V, Y, _, Coverage),
#     not(decision(C, V, Y, acceptReported, _, _, _)),
#     Coverage > 100,
#     Reason = 100,
#     concat_atom(['Reported data excluded because ', Coverage,
#     ' percent greater than 100 percent. '], Expl).

# index = Rep.Cov > 100
# accept = Decisions[Decisions$Dec == "acceptReported", ]
# index[accept$Y, accept$V] = FALSE
# Expl[which(index)] = sprintf(
#     "%sReported data excluded because %i percent greater than 100 percent. ", 
#     Expl[which(index)], Rep.Cov[which(index)])

# reported_reason_to_exclude(C, V, Y, _, Expl) :-
#     reported(C, V, Y, _, Coverage),
#     not(decision(C, V, Y, acceptReported, _, _, _)),
#     Prec is Y - 1,
#     Succ is Y + 1,
#     reported(C, V, Prec, _, CoveragePrec),
#     reported(C, V, Succ, _, CoverageSucc),
#     sawtooth_threshold(Threshold),
#     Coverage - CoveragePrec > Threshold,
#     Coverage - CoverageSucc > Threshold,
#     concat_atom(
#       [ 'Reported data excluded due to an increase from ', CoveragePrec,
#         ' percent to ', Coverage, ' percent with decrease ', CoverageSucc,
#         ' percent. '
#       ], Expl).
#
# reported_reason_to_exclude(C, V, Y, _, Expl) :-
#     reported(C, V, Y, _, Coverage),
#     not(decision(C, V, Y, acceptReported, _, _, _)),
#     Prec is Y - 1,
#     Succ is Y + 1,
#     reported(C, V, Prec, _, CoveragePrec),
#     reported(C, V, Succ, _, CoverageSucc),
#     sawtooth_threshold(Threshold),
#     CoveragePrec - Coverage > Threshold,
#     CoverageSucc - Coverage > Threshold,
#     concat_atom(
#       [ 'Reported data excluded due to decline in reported coverage from ',
#         CoveragePrec, ' percent to ', Coverage, ' percent with increase to ',
#         CoverageSucc,' percent. '], Expl).

index = which(Rej.Info != "")
Expl[index] = sprintf("%s%s", Expl[index], Rej.Info[index])

# explanation(C, V, Y, Expl) :-
#     decision(C, V, Y, comment, Expl, _, _).

index = Decisions[Decisions$Dec == "comment", ]
if(nrow(index))
  for(i in 1:nrow(index))
    Expl[index$Y[i], index$V[i]] =
      sprintf("%s%s", Expl[index$Y[i], index$V[i]], index$Info[i])

# explanation(C, V, Y, Expl) :-
#     decision(C, V, Y, acceptSurvey, Expl, _, _).

index = Decisions[Decisions$Dec == "acceptSurvey", ]
if(nrow(index))
  for(i in 1:nrow(index))
    Expl[index$Y[i], index$V[i]] =
      sprintf("%s%s", Expl[index$Y[i], index$V[i]], index$Info[i])

# explanation(C, V, Y, Expl) :-
#     decision(C, V, Y, acceptReported, Expl, _, _).

index = Decisions[Decisions$Dec == "acceptReported", ]
if(nrow(index))
  for(i in 1:nrow(index))
    Expl[index$Y[i], index$V[i]] =
      sprintf("%s%s", Expl[index$Y[i], index$V[i]], index$Info[i])

# explanation(C, V, Y, Expl) :-
#     decision(C, V, Y, ignoreGov, Expl, _, _).

index = Decisions[Decisions$Dec == "ignoreGov", ]
if(nrow(index))
  for(i in 1:nrow(index))
    Expl[index$Y[i], index$V[i]] =
      sprintf("%s%s", Expl[index$Y[i], index$V[i]], index$Info[i])

Text = YV_char
Text[] = sprintf("%s %s %s %s", Info, Expl, Change, GoC.Expl)

Vaccine = Vn
Year = Yn
VY = expand.grid(Year, Vaccine, stringsAsFactors=FALSE)

include = Ereq & !is.na(Bounded)
VY = cbind(Y=VY$Var1[include], V=VY$Var2[include])

Table = data.frame(
    Country=Country,
    ProductionDate=Date,
    ISOCountryCode=Code,
    Vaccine=VY[, "V"],
    Year=VY[, "Y"],
    WUENIC=Bounded[VY],
    WUENICPreviousRevision=Legacy[VY],
    GradeOfConfidence=GoC[VY],
    AdministrativeCoverage=Admin[VY],
    GovernmentEstimate=Gov[VY],
    ReportedCoverage=Rep.Cov[VY],
    ChildrenVaccinated=Vaccinated[VY],
    ChildrenInTarget=Target[VY],
    BirthsUNPD=Births[VY[, "Y"]],
    SurvivingInfantsUNPD=Surviving[VY[, "Y"]],
    ReportedTimeSeries=TS.Cov[VY],
    ReportedTimeSeriesSource=TS.Src[VY],
    SurveyInformation=Svy.Cov[VY],
    Rule=Rule[VY],
    Comment=Text[VY])

dir.create("out", showWarnings=FALSE, recursive=TRUE)
write.table(Table, sprintf("out/%s.txt", ccode), quote=FALSE, row.names=FALSE, sep="\t", na="")

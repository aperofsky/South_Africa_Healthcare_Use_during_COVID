###############################################################
## South Africa Netcare Inpatient Excess Respiratory Rate -National
## Amanda Perofsky amanda.perofsky@nih.gov
###############################################################
## note: code will not run because the input data require DSAs with
## NICD and Netcare.


# rm(list=ls())

list.of.packages <- c("dplyr","data.table","tidyverse","padr","ggplot2","cowplot",
                      "forecast","ciTools","lubridate","cdcfluview","zoo","sweep","tidyquant","magrittr","viridis",
                      "ggExtra","rcartocolor","wesanderson","RColorBrewer")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

# set working directory to South_Africa_ILI folder
setwd("~/South_Africa_ILI/")
## Output will be in "Manuscript_CSVs" folder within "South_Africa_ILI"
## Figures will be in "Manuscript_Excess_Cases_Figures" folder within "South_Africa_ILI"

#first week: 2020-03-01
last_date = "2021-09-26" #specify last epidemic week here
mmwr_week(as.Date(last_date)) # 2021 week 39
phases = data.frame(epi_dates = seq(from=as.Date("2020-03-01"),to=as.Date(last_date),by="week"),
                    lockdown=c("pre-lockdown","pre-lockdown","pre-lockdown","pre-lockdown",
                               "phase 5","phase 5","phase 5","phase 5","phase 5",
                               "phase 4","phase 4","phase 4","phase 4",
                               "phase 3","phase 3","phase 3","phase 3",
                               "phase 3","phase 3","phase 3","phase 3",
                               "phase 3","phase 3","phase 3","phase 2","phase 2",
                               "phase 2","phase 2","phase 2","phase 1","phase 1",
                               "phase 1","phase 1","phase 1","phase 1","phase 1",
                               "phase 1","phase 1","phase 1","phase 1","phase 1",
                               "phase 1","phase 1","adj. phase 3","adj. phase 3","adj. phase 3",
                               "adj. phase 3","adj. phase 3","adj. phase 3","adj. phase 3","adj. phase 3",
                               "adj. phase 3","adj. phase 1","adj. phase 1","adj. phase 1","adj. phase 1",
                               "adj. phase 1","adj. phase 1","adj. phase 1","adj. phase 1","adj. phase 1",
                               "adj. phase 1","adj. phase 1","adj. phase 1","adj. phase 1","adj. phase 2","adj. phase 2",
                               "adj. phase 3.2", "adj. phase 3.2","adj. phase 4", "adj. phase 4",
                               "adj. phase 4","adj. phase 4","adj. phase 3.3","adj. phase 3.3",
                               "adj. phase 3.3","adj. phase 3.3","adj. phase 3.3","adj. phase 3.3",
                               "adj. phase 3.3","adj. phase 2.2","adj. phase 2.2","adj. phase 2.2"))

phases

## national covid cases
covid_cases = fread("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
sa_covid_cases = covid_cases %>% filter(`Country/Region`=="South Africa")
sa_covid_cases = sa_covid_cases %>% pivot_longer(names_to="date",cols=c(`1/22/20`:last(colnames(sa_covid_cases))),values_to="confirmed")
sa_covid_cases$date = as.Date(sa_covid_cases$date,format = "%m/%d/%y")
sa_covid_cases = bind_cols(sa_covid_cases,cdcfluview::mmwr_week(sa_covid_cases$date)[1:2])
sa_covid_cases$cdc_date = cdcfluview::mmwr_week_to_date(sa_covid_cases$mmwr_year,sa_covid_cases$mmwr_week)

sa_covid_cases = 
  sa_covid_cases %>%
  mutate(new_confirmed=confirmed-dplyr::lag(confirmed,lag=1,order_by = date))%>%
  mutate(new_confirmed=replace_na(new_confirmed,0))%>%
  dplyr::select(date,mmwr_year,mmwr_week,cdc_date,new_confirmed)

weekly_cum = sa_covid_cases %>% 
  group_by(mmwr_year,mmwr_week,cdc_date) %>%
  summarize(weekly_cases = sum(new_confirmed,na.rm = T)) %>%
  ungroup()
names(weekly_cum)[1:2]<- c("year","week")

###############################################################
## three distinct surveillance systems: SARI, ILI Public, ILI Watch
## aggegate RSV+ and flu+ across all systems
###############################################################
## import most recent data draw down (weekly ILI and SARI cases)
dir<-"data"
folder = list.files(path = dir, pattern = "SA_Data",full.names = TRUE) %>%
  magrittr::extract(which.max(file.mtime(.)))
folder

dl_date = as.Date(strsplit(folder,"data/SA_Data_")[[1]][2],"%m_%d_%Y")
epi_week = mmwr_week(dl_date)
epi_date = mmwr_week_to_date(epi_week$mmwr_year,epi_week$mmwr_week)

sari = read.csv(paste0(folder,"/SARI - 2015-2021.csv"))
sari$dataset = "SARI"
sari$weekdate = parse_date_time(sari$weekdate, orders = c('dmy'))
sari$weekdate = as.Date(sari$weekdate)
sari = sari %>% mutate(weekdate=if_else(is.na(weekdate),as.Date(mmwr_week_to_date(year,week)),as.Date(weekdate)))

ilipub = read.csv(paste0(folder,"/ILI - Public Clinics - 2015-2021.csv"))
ilipub$dataset = "ILI Public"
ilipub$weekdate = parse_date_time(ilipub$weekdate, orders = c('dmy'))
ilipub$weekdate = as.Date(ilipub$weekdate)
ilipub = ilipub %>% mutate(weekdate=if_else(is.na(weekdate),as.Date(mmwr_week_to_date(year,week)),as.Date(weekdate)))

ili = bind_rows(ilipub,sari) 

ili.sa=ili %>%
  group_by(weekdate)  %>%
  summarize(sum(flu,na.rm = T), 
            sum(rsv,na.rm = T), 
            sum(cov,na.rm = T),
            sum(flursvcovneg,na.rm = T), 
            sum(total,na.rm = T))
ili.sa=data.table(ili.sa)
names(ili.sa)=c("weekdate","flu","rsv","cov", "neg","total")

ili.sa = 
  ili.sa %>%
  padr::pad()%>%
  mutate(flu = na.approx(flu,maxgap = 6,na.rm = F),
         rsv = na.approx(rsv,maxgap = 6,na.rm = F),
         total = na.approx(total,maxgap = 6,na.rm = F)) %>%
  mutate(flu = ifelse(is.na(flu),0,flu),
         rsv = ifelse(is.na(rsv),0,rsv),
         total = ifelse(is.na(total),0,total)) %>%
  mutate(per_rsv = rsv/total,
         per_flu = flu/total)

#one week lags of per_flu and per_rsv
ili.sa.lag = ili.sa %>%
  mutate(one_week_lag_per_flu = lag(per_flu,order_by = weekdate,n = 1),
         one_week_lag_per_rsv = lag(per_rsv,order_by = weekdate,n = 1)) %>%
  fill(one_week_lag_per_flu:one_week_lag_per_rsv,.direction = "up") #fill in first week of 2015 using the second week's value

ili.sa.lag = ili.sa.lag %>%
  fill(.direction = "down")

ili.sa.lag$year = mmwr_week(ili.sa.lag$weekdate)[,1]
ili.sa.lag$week = mmwr_week(ili.sa.lag$weekdate)[,2]

prev_season = ili.sa.lag %>%
  filter(weekdate>=as.Date("2018-12-28") & weekdate<=as.Date("2019-12-28"))%>%
  dplyr::select(week,per_flu,per_rsv)
names(prev_season)[2:3]<- c("flu_2019","rsv_2019")
head(prev_season)
class(prev_season$week)
class(ili.sa.lag$week)
nrow(prev_season)

ili.sa.lag=left_join(ili.sa.lag,prev_season,by="week")
# replaces data from 2020 onwards with data from 2019
ili.sa.lag = ili.sa.lag %>% 
  mutate(avg_per_flu = ifelse(weekdate>=as.Date("2020-03-01"),flu_2019,per_flu),
         avg_per_rsv = ifelse(weekdate>=as.Date("2020-03-01"),rsv_2019,per_rsv))%>%
  ungroup()

ili.sa.lag=
  ili.sa.lag %>%
  mutate(avg_per_flu = rollmean(x = avg_per_flu, 4, align = "center", fill = NA),
         avg_per_rsv = rollmean(x = avg_per_rsv, 4, align = "center", fill = NA))%>%
  fill(c(avg_per_flu,avg_per_rsv),.direction = "updown")
#################################
## All cause respiratory admissions
#################################
## edit code to most recent data drawdown (publication data does not include any nowcasted admissions)
backfill <-"National_Backfill_Analysis"
subfile = list.files(path = backfill, pattern = "netcare_inpt_nowcast_df_as_of_",full.names = TRUE) 
national_file = subfile %>% extract(which.max(file.mtime(.)))
load(national_file)
netcare_inpt.sa = netcare_inpt_forecasting_df

netcare_inpt.sa = netcare_inpt.sa %>%
  mutate(total_res = if_else(is.na(total_res),0,total_res),
         total_covid = if_else(is.na(total_covid),0,total_covid),
         total_visits = if_else(is.na(total_visits),0,total_visits))%>%
  mutate(total_res_cov = total_res+total_covid)%>%
  mutate(res_cov_rate = total_res_cov/total_visits)

############################################################################
## combine admissions data with lab-confirmed pathogen data
############################################################################

netcare_inpt.sa= left_join(netcare_inpt.sa,ili.sa.lag %>% 
                             dplyr::select(weekdate,flu,rsv,neg,total,
                                           per_flu,per_rsv,avg_per_flu,
                                           avg_per_rsv,one_week_lag_per_flu,
                                           one_week_lag_per_rsv),
                           by=c("weekdate"))
netcare_inpt.sa= tidyr::complete(netcare_inpt.sa,weekdate)

netcare_inpt.sa = netcare_inpt.sa %>%
  fill(c(flu,rsv,neg,total,per_flu,per_rsv,avg_per_flu,
         avg_per_rsv,one_week_lag_per_flu,one_week_lag_per_rsv), .direction = "down")%>%
  ungroup()


netcare_inpt.sa.covid = left_join(netcare_inpt.sa,weekly_cum %>% 
                                    dplyr::select(-cdc_date),by=c("year","week")) %>%
                        mutate(weekly_cases = ifelse(is.na(weekly_cases),0,weekly_cases))

netcare_inpt.sa.covid = netcare_inpt.sa.covid %>% filter(weekdate<=last_date)
###################################################################
## national all-cause respiratory admissions time series models - counts
##################################################################

agegrp.model.netcare.in.counts <- function(case_factor = 20,outcome="total_res_cov",k1=6,k2=9, label="Respiratory Encounters",plot.title="Inpatient") {  
  
  df = netcare_inpt.sa.covid %>% 
    filter(weekdate>=as.Date("2016-01-01")) %>%
    droplevels()
  max_date = max(df$weekdate)
  
  df$one_week_lag_per_flu = scale(df$one_week_lag_per_flu)[,1]
  df$one_week_lag_per_rsv = scale(df$one_week_lag_per_rsv)[,1]
  
  # max(df$weekdate)
  # min(df$weekdate)
  # nrow(df)
  
  test_dates = df %>% filter(weekdate>as.Date("2020-02-23")) %>% distinct(weekdate)
  covid_weekly = df %>% filter(weekdate>as.Date("2020-02-23")) %>% dplyr::select(weekly_cases) 
  covid_netcare = df %>% filter(weekdate>as.Date("2020-02-23")) %>% dplyr::select(total_covid) 
  
  df_ts = msts(df,seasonal.periods = c(52.18,26.09),start=2016)
  tail(df_ts)
  
  res_rate <- df_ts[,outcome] #entire time series
  
  df_ts_train = window(df_ts,end=c(2020,9)) #training set
  df_ts_test = window(df_ts,start=c(2020,10)) #test set
  
  test_dates =df %>% filter(weekdate>as.Date("2020-02-23"))

  
  res_rate_train = df_ts_train[,outcome]
  res_rate_test = df_ts_test[,outcome]
  
  # rsv and flu positive samples
  lab_train = df_ts_train[,c("avg_per_flu","avg_per_rsv")]
  z_train =  fourier(df_ts_train, K=c(k1,k2))
  xreg_train = as.matrix(data.frame(z_train,lab_train))
  colnames(xreg_train)

  lab_test = df_ts_test[,c("avg_per_flu","avg_per_rsv")]
  z_test =  fourier(df_ts_test, K=c(k1,k2))
  xreg_test = as.matrix(data.frame(z_test,lab_test))

  # # flu+, rsv+ current week
  fit_res_rate_flu_rsv_current <- auto.arima(res_rate_train,xreg=xreg_train,seasonal=F, lambda=0,stepwise = F,approximation = F,parallel = T)

  test_len = length(res_rate_test) 
  
  forecast_res_rate = fit_res_rate_flu_rsv_current  %>% forecast(h=test_len,xreg=xreg_test,bootstrap=T,npaths=1000)

  vec_length = length(head(res_rate,-test_len))
  ts_length = length(res_rate)
  
  
  date_seq<-format(lubridate::date_decimal(as.numeric(row.names(as.data.frame(forecast_res_rate)))),"%Y-%m-%d")
  test_date_vec = as.numeric(time(res_rate_test))
  train_date_vec = as.numeric(time(res_rate_train))
  length(test_date_vec)
  
  covid_df = data.frame(index = date_seq,key="COVID cases",
                        value=covid_netcare * case_factor)
  test_df =data.frame(index=test_date_vec,cases=as.vector(res_rate_test))
  ts_df = data.frame(index=as.numeric(time(df_ts)),cases=df_ts[,outcome])
  
  unexplained_cases = as.vector(res_rate_train - fitted(fit_res_rate_flu_rsv_current))
  unexplained_cases[unexplained_cases < 0] <- 0
  unexplained_cases_forecast = res_rate_test - as.data.frame(forecast_res_rate)[,1]
  unexplained_cases_forecast[unexplained_cases_forecast < 0] <- 0
  unexplained_cases_combined = c(unexplained_cases,unexplained_cases_forecast)
  predicted = c(fitted(fit_res_rate_flu_rsv_current),as.data.frame(forecast_res_rate)[,1])
  
  averted_cases = as.vector(res_rate_train - fitted(fit_res_rate_flu_rsv_current))
  averted_cases[averted_cases > 0] <- 0
  averted_cases = abs(averted_cases)
  averted_cases_binary = averted_cases
  averted_cases_binary[averted_cases_binary>0]<- 1
  
  averted_cases_forecast = res_rate_test - as.data.frame(forecast_res_rate)[,1]
  averted_cases_forecast[averted_cases_forecast > 0] <- 0
  averted_cases_binary_forecast = averted_cases_forecast
  averted_cases_binary_forecast[averted_cases_binary_forecast<0]<- 1
  averted_cases_combined = c(averted_cases_binary,averted_cases_binary_forecast)
  
  
  
  plot.df = data.frame(covid = df$total_covid,
                       unexplained = unexplained_cases_combined,
                       obs = df$total_res_cov,pi=df$total_pi,
                       pred = predicted,
                       date = df$weekdate,
                       averted_cases = averted_cases_combined)
  plot.df =  plot.df %>%
    mutate(cases_not_due_to_covid = obs-covid)
  
  plot.df2 = plot.df %>%
    mutate(prop.covid = covid/obs,
           prop.pneu_and_influ = pi/obs,
           prop.resp_not_PI_or_covid = 1 - ((pi+covid)/obs)) %>%
    mutate(excess.covid = unexplained*prop.covid,
           excess.pneu_and_influ = unexplained*prop.pneu_and_influ,
           excess.non_viral.resp = unexplained*prop.resp_not_PI_or_covid) %>%
    mutate(excess.covid = ifelse(excess.covid<0,0,excess.covid),
           excess.pneu_and_influ = ifelse(excess.pneu_and_influ<0,0,excess.pneu_and_influ),
           excess.non_viral.resp = ifelse(excess.non_viral.resp<0,0,excess.non_viral.resp))
  
  tidy_forecast = sw_sweep(forecast_res_rate, fitted = TRUE)
  forecast_df = tidy_forecast %>% filter(key=="forecast")
  nrow(forecast_df)
  length(test_dates$weekdate)
  forecast_df = data.frame(forecast_df,date=test_dates$weekdate)
  forecast_df$date = as.Date(forecast_df$date)
  
  plot.df3 = left_join(plot.df2,forecast_df %>% dplyr::select(-index,-key,-value),by="date")
  tail(plot.df3)
  
  plot.df4 = plot.df3 %>%
    mutate(max_obs = pmax(obs,hi.95),
           min_obs = pmin(obs,lo.95))
  
  stack_plot_func <- function(ds, legend=T, scale.plot=1,legend.loc='topleft', set.ymax=0){
    ds$band1 <- ds$pred +  ds$excess.covid
    ds$band1.full <- ds$pred +  ds$excess.covid
    
    ds$band2 <- ds$band1 +  ds$excess.non_viral.resp
    
    ds$band2[ds$band2> ds$obs] <- ds$obs[ds$band2> ds$obs]
    
    ds <- ds[!is.na(ds$pred),]
    initial.base <- ds$pred[1]
    
    #yrange.pneu <- range(c(ds$pred, ds$obs, ds$band2),0)
    if(set.ymax==0){
      yrange.pneu <- c(0, initial.base*scale.plot)
    }else{
      yrange.pneu <- c(0, set.ymax)
    }
    
    
    nb.cols <- 11
    cols <- brewer.pal(n=nb.cols,"RdYlBu")
    
    ## data are weekly so dates correspond to the epidemic weeks of each alert level
    p = ggplot(ds %>% filter(date>as.Date("2018-01-01")))+
      ylim(yrange.pneu)+
      ylab(label)+
      xlab("Date")+
      scale_x_date(date_breaks = "1 year",date_labels = "%Y",expand=c(0,0))+
      annotate("rect",xmin= as.Date("2020-03-01"),
               xmax =as.Date("2020-03-29"),
               ymin = -Inf,
               ymax = Inf, fill = "light blue", alpha = 0.4)+
      annotate("rect",xmin= as.Date("2020-03-29"),
               xmax =as.Date("2020-05-03"),
               ymin = -Inf,
               ymax = Inf, fill = cols[1], alpha = 0.4)+##phase 5
      annotate("rect",xmin= as.Date("2020-05-03"),
               xmax = as.Date("2020-05-30"),
               ymin = -Inf,
               ymax = Inf, fill = cols[2], alpha = 0.4)+#phase 4
      annotate("rect",xmin= as.Date("2020-05-30"),
               xmax = as.Date("2020-08-16"),
               ymin = -Inf,
               ymax = Inf, fill = cols[3], alpha = 0.4)+#phase 3
      annotate("rect",xmin= as.Date("2020-08-16"),
               xmax =as.Date("2020-09-20"),
               ymin = -Inf,
               ymax = Inf, fill = cols[4], alpha = 0.4) +#phase 2
      annotate("rect",xmin= as.Date("2020-09-20"),
               xmax = as.Date("2020-12-27"),
               ymin = -Inf,
               ymax = Inf, fill = cols[5], alpha = 0.4) +#phase 1
      annotate("rect",xmin= as.Date("2020-12-27"),
               xmax = as.Date("2021-02-28"),
               ymin = -Inf,
               ymax = Inf, fill = cols[3], alpha = 0.4) +#adj. phase 3
      annotate("rect",xmin= as.Date("2021-02-28"),
               xmax = as.Date("2021-05-30"),
               ymin = -Inf,
               ymax = Inf, fill = cols[5], alpha = 0.4) +#adj. phase 1
      annotate("rect",xmin= as.Date("2021-05-30"),
               xmax = as.Date("2021-06-13"),
               ymin = -Inf,
               ymax = Inf, fill = cols[4], alpha = 0.4) +#adj. phase 2
      annotate("rect",xmin= as.Date("2021-06-13"),
               xmax = as.Date("2021-06-27"),
               ymin = -Inf,
               ymax = Inf, fill = cols[3], alpha = 0.4) +#adj. phase 3
      annotate("rect",xmin= as.Date("2021-06-27"),
               xmax = as.Date("2021-07-25"),
               ymin = -Inf,
               ymax = Inf, fill = cols[2], alpha = 0.4) +#adj. phase 4
      annotate("rect",xmin= as.Date("2021-07-25"),
               xmax = as.Date("2021-09-12"),
               ymin = -Inf,
               ymax = Inf, fill = cols[3], alpha = 0.4) +#adj. phase 3     
      annotate("rect",xmin= as.Date("2021-09-12"),
               xmax = max(ds$date),
               ymin = -Inf,
               ymax = Inf, fill = cols[4], alpha = 0.4) +#adj. phase 2     
      geom_ribbon(data = ds %>% filter(date>=as.Date("2020-03-01")),aes(x=date,ymin=pred,ymax=max_obs,fill='Excess cases'), alpha = 1)+
      geom_ribbon(data=ds %>% filter(date>=as.Date("2020-03-01")),aes(x=date,ymin=obs,ymax=max_obs,fill='Averted cases'))+
      geom_ribbon(data=ds %>% filter(date>as.Date("2018-01-01") & date<=as.Date("2020-03-01")), aes(x=date,ymin=0,ymax=pred),fill="#EAE7E5")+
      geom_ribbon(data=forecast_df,aes(x=date,ymin = lo.95, ymax = hi.95),
                  fill = "#D5DBFF", color = NA, size = 0,alpha=0.8) +
      geom_ribbon(data=forecast_df,aes(x=date,ymin = lo.80, ymax = hi.80, fill = key),
                  fill = "#596DD5", color = NA, size = 0, alpha = 0.5) +
      geom_ribbon(data = ds %>% filter(date >= as.Date("2020-03-01")), aes(x=date,ymin=0,ymax=min_obs),fill="#EAE7E5")+
      geom_line(aes(x=date,y=obs,color="Observed respiratory\nadmissions"),lwd=1.5)+
      geom_line(aes(x=date,y=covid,color="COVID coded (U07)\nadmissions"),lty="dashed",lwd=1)+
      geom_line(aes(x=date,y=pred,color="Seasonal\nbaseline"),lwd=1.5)+
      geom_vline(xintercept = as.Date("2020-03-01"),lty="dashed",color="black",lwd=1)+
      scale_fill_manual(name="",
                        values= c(wes_palette("Moonrise3")[3],wes_palette("Zissou1")[3]),
                         labels= c('Averted cases','Excess cases'))+
      scale_color_manual(name="", values=c("#F8766D","black","#596DD5"))+
      theme_bw(base_size = 12)+
      background_grid(major="none",minor="none")+
      scale_y_continuous(expand=c(0,Inf))+
      labs(title = plot.title) +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.text.align=0.5,
            axis.text.y = element_text(size = 10),
            legend.position = "bottom",
            legend.box="horizontal", 
            legend.margin=margin(),
            text = element_text(family = "Arial"))
    print(p)
  }
  p = stack_plot_func(ds=plot.df4,legend = T,set.ymax = max(plot.df2$obs)+50) 
  
  excess_upper = as.vector(res_rate_test[1:test_len]-forecast_res_rate$upper[1:test_len,2])
  excess_lower = as.vector(res_rate_test[1:test_len]-forecast_res_rate$lower[1:test_len,2])
  excess_point = as.vector(res_rate_test[1:test_len]-forecast_res_rate$mean[1:test_len])
  
  excess_point_percent_diff = 100*(as.vector(res_rate_test[1:test_len]-forecast_res_rate$mean[1:test_len])/forecast_res_rate$mean[1:test_len])
  excess_lower_percent_diff = 100*(as.vector(res_rate_test[1:test_len]-forecast_res_rate$lower[1:test_len,2])/forecast_res_rate$lower[1:test_len,2])
  excess_upper_percent_diff = 100*(as.vector(res_rate_test[1:test_len]-forecast_res_rate$upper[1:test_len,2])/forecast_res_rate$upper[1:test_len,2])
  
  df = data.frame(encounter_type="inpatient",
                  dates=test_dates$weekdate,
                  excess_point = round(excess_point),
                  excess_ulim=round(excess_upper),
                  excess_llim = round(excess_lower),
                  excess_point_per_diff = excess_point_percent_diff,
                  excess_ulim_per_diff = excess_upper_percent_diff,
                  excess_llim_per_diff=excess_lower_percent_diff,
                  covid_admissions = covid_netcare,
                  national_covid_cases = covid_weekly,
                  pred_baseline = as.vector(forecast_res_rate$mean[1:test_len]),
                  obs_cases = as.vector(res_rate_test[1:test_len]))
  df
  nrow(df)
  
  df = bind_cols(df,phases)
  
  return(list(df, p))
}

load("Manuscript_ARIMA/netcare_inpt_national_fourier_terms.Rdata")
inpt_national_fourier

p = agegrp.model.netcare.in.counts(k1=inpt_national_fourier$k1,k2=inpt_national_fourier$k2, outcome="total_res_cov",label="Admissions",plot.title="Inpatient")

max_date = max(netcare_inpt.sa.covid$weekdate)

legend_fill <- get_legend(
  p[[2]] +
    guides(fill= guide_legend(nrow = 1),
           lty=guide_legend(nrow=1)))


resp_visit_counts_fig_arima_netcare_in = p[[2]]+theme(legend.position = "none")
save_plot(filename = paste0("Manuscript_Excess_Cases_Figures/netcare_inpt_arima_flu_resp_counts_national_as_of_",max_date,".png"),resp_visit_counts_fig_arima_netcare_in, base_width = 10, base_height = 5,dpi=600)


per_diff_resp = p[[1]]
write.csv(per_diff_resp,file = paste0("Manuscript_CSVs/national_netcare_inpt_all_cause_respiratory_excess_estimates_arima_as_of_",max_date,".csv"),row.names = F)


nb.cols <- 11
cols <- brewer.pal(n=nb.cols,"RdYlBu")


per_diff_plot = ggplot(per_diff_resp)+
  geom_hline(aes(yintercept=0),lty="dashed",color="black")+
  xlab("Date")+
  ylab("Percent Difference from Baseline")+
  theme_bw(base_size = 12)+
  geom_vline(xintercept = as.Date("2020-03-01"),lty="dashed",color="black",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2020-03-29")),color=cols[1],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2020-05-03")),color=cols[2],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2020-05-30")),color=cols[3],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2020-08-16")),color=cols[4],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2020-09-20")),color=cols[5],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2020-12-27")),color=cols[3],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2021-02-28")),color=cols[5],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2021-05-30")),color=cols[4],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2021-06-13")),color=cols[3],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2021-06-27")),color=cols[2],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2021-07-25")),color=cols[3],lty="dashed",lwd=1)+
  geom_vline(aes(xintercept=as.Date("2021-09-12")),color=cols[4],lty="dashed",lwd=1)+
  annotate("rect",xmin= as.Date("2020-03-01"),
           xmax =as.Date("2020-03-29"),
           ymin = -Inf,
           ymax = Inf, fill = "light blue", alpha = 0.4)+
  annotate("rect",xmin= as.Date("2020-03-29"),
           xmax =as.Date("2020-05-03"),
           ymin = -Inf,
           ymax = Inf, fill = cols[1], alpha = 0.4)+##phase 5
  annotate("rect",xmin= as.Date("2020-05-03"),
           xmax = as.Date("2020-05-30"),
           ymin = -Inf,
           ymax = Inf, fill = cols[2], alpha = 0.4)+#phase 4
  annotate("rect",xmin= as.Date("2020-05-30"),
           xmax = as.Date("2020-08-16"),
           ymin = -Inf,
           ymax = Inf, fill = cols[3], alpha = 0.4)+#phase 3
  annotate("rect",xmin= as.Date("2020-08-16"),
           xmax =as.Date("2020-09-20"),
           ymin = -Inf,
           ymax = Inf, fill = cols[4], alpha = 0.4) +#phase 2
  annotate("rect",xmin= as.Date("2020-09-20"),
           xmax = as.Date("2020-12-27"),
           ymin = -Inf,
           ymax = Inf, fill = cols[5], alpha = 0.4) +#phase 1
  annotate("rect",xmin= as.Date("2020-12-27"),
           xmax = as.Date("2021-02-28"),
           ymin = -Inf,
           ymax = Inf, fill = cols[3], alpha = 0.4) +#adj. phase 3
  annotate("rect",xmin= as.Date("2021-02-28"),
           xmax = as.Date("2021-05-30"),
           ymin = -Inf,
           ymax = Inf, fill = cols[5], alpha = 0.4) +#adj. phase 1
  annotate("rect",xmin= as.Date("2021-05-30"),
           xmax = as.Date("2021-06-13"),
           ymin = -Inf,
           ymax = Inf, fill = cols[4], alpha = 0.4) +#adj. phase 2
  annotate("rect",xmin= as.Date("2021-06-13"),
           xmax = as.Date("2021-06-27"),
           ymin = -Inf,
           ymax = Inf, fill = cols[3], alpha = 0.4) +#adj. phase 3
  annotate("rect",xmin= as.Date("2021-06-27"),
           xmax = as.Date("2021-07-25"),
           ymin = -Inf,
           ymax = Inf, fill = cols[2], alpha = 0.4) +#adj. phase 4
  annotate("rect",xmin= as.Date("2021-07-25"),
           xmax = as.Date("2021-09-12"),
           ymin = -Inf,
           ymax = Inf, fill = cols[3], alpha = 0.4) +#adj. phase 3
  annotate("rect",xmin= as.Date("2021-09-12"),
           xmax = max(per_diff_resp$date),
           ymin = -Inf,
           ymax = Inf, fill = cols[4], alpha = 0.4) +#adj. phase 2
  geom_point(aes(x=dates,y=excess_point_per_diff))+
  geom_errorbar(aes(x=dates,ymin=excess_llim_per_diff, ymax=excess_ulim_per_diff), colour="black", width=.1)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(family = "Arial"))+
  scale_x_date(date_breaks = "3 months",date_labels = "%b %y",expand=c(0.02,0.02))
per_diff_plot

resp_figure = plot_grid(p[[2]],per_diff_plot,nrow=2,rel_heights = c(1,0.75))
resp_figure
save_plot(resp_figure,filename = paste0("Manuscript_Excess_Cases_Figures/national_resp_netcare_inpatient_excess_as_of_",max_date,".png"),base_width = 12, base_height = 12,dpi=600)


require 'rbbt-util'
require 'rbbt/workflow'

module NodeFeatures
  extend Workflow

  
  input :target_feature, :string, "name of feature to predict", 'anyCells'
  input :feature_files, :array, "Feature files to load", nil, :nofile => true
  input :balance_data, :booles, "Make balanced class priors", false
  task :importance => :tsv do |target_feature,feature_files, balance_data|
    data = Rbbt.data["Classification_influential_nodes.txt"].tsv(:fields => [target_feature], :type => :list)
    feature_files = Rbbt.data.glob("*.txt").collect{|f| f.find} if feature_files.nil?
    feature_files = feature_files - [Rbbt.data["Classification_influential_nodes.txt"].find]

    feature_files.each do |file|
      ftsv = TSV.open(file)
      ftsv.fields = ftsv.fields.collect{|f| f.sub(":",'_').strip}
      ffields = ftsv.fields - ["HGNC.symbols"]
      double = Misc.counts(ffields).select{|f,c| c > 1}.collect{|f,c| f}
      ffields = ffields - double - data.fields 

      ftsv.process "Drug.target" do |value|
        value = value.first if Array === value
        value == "True" ? [1] : [0]
      end if ftsv.fields.include? "Drug.target"

      data = data.attach ftsv, :fields => ffields
    end

    factor_fields = []
    data.fields.each do |field|
      next if field == target_feature
      if data.column(field).values.flatten.compact.uniq.sort == ["0", "1"]
        factor_fields << field
        data.process field do |v|
          v == "0" ? 'not-present' : 'present'
        end
      end
    end 

    require 'rbbt/util/R'
    Open.mkdir files_dir
    Open.write(file('factor_fields'), factor_fields * "\n")
    Open.write(file('data'), data.to_s)
    data.R <<-EOF
rbbt.require('randomForest')
names(data) = make.names(names(data))
factor.fields = scan('#{file('factor_fields')}', what=character(), sep="\n")
factor.fields = make.names(factor.fields)

for (f in factor.fields){
  data[,f] = as.factor(data[,f])
}

data[is.na(data)] = -1
data = data[data$#{target_feature} != -1,]
if (#{balance_data.to_s.upcase}){
  m <- randomForest(as.factor(#{target_feature}) ~ ., data, na.action=na.omit, classwt=c(0.5,0.5))
}else{
  m <- randomForest(as.factor(#{target_feature}) ~ ., data, na.action=na.omit)
}
save(m, file = '#{self.file('model.Rdata')}')
data = NULL
    EOF

    confusion = data.R <<-EOF
rbbt.require('randomForest')
load(file = '#{self.file('model.Rdata')}')
data = m$confusion
    EOF

    Open.write(file('confusion.tsv'), confusion.to_s)

    importance = data.R <<-EOF
rbbt.require('randomForest')
load(file = '#{self.file('model.Rdata')}')
data = m$importance
    EOF

    Open.write(file('importance.tsv'), importance.to_s)

    err = data.R <<-EOF
rbbt.require('randomForest')
load(file = '#{self.file('model.Rdata')}')
data = m$err.rate
    EOF

    log :err_rate, err[err.keys.last].first
    set_info :err_rate, err[err.keys.last].first
    Open.write(file('err_rate.tsv'), err.to_s)

    confusion
  end

  input :target_feature, :string, "name of feature to predict", 'anyCells'
  input :feature_files, :array, "Feature files to load", nil, :nofile => true
  task :cforest_importance => :tsv do |target_feature,feature_files|
    data = Rbbt.data["Classification_influential_nodes.txt"].tsv(:fields => [target_feature], :type => :list)
    feature_files = Rbbt.data.glob("*.txt").collect{|f| f.find} if feature_files.nil?
    feature_files = feature_files - [Rbbt.data["Classification_influential_nodes.txt"].find]

    feature_files.each do |file|
      ftsv = TSV.open(file)
      ffields = ftsv.fields - ["HGNC.symbols"]
      ffields.reject!{|f| f != f.strip}
      data = data.attach ftsv, :fields => ffields
    end


    data.fields.each do |field|
      next if field == target_feature
      if data.column(field).values.flatten.compact.uniq.sort == ["0", "1"]
        data.process field do |v|
          v == "0" ? 'not-present' : 'present'
        end
      end
    end


    require 'rbbt/util/R'
    Open.mkdir files_dir

    Open.write(file('data'), data.to_s)

    data.R <<-EOF
rbbt.require('party')
names(data) = make.names(names(data))

data[is.na(data)] = -1
data = data[data$#{target_feature} != -1,]
data = as.data.frame(lapply(data, as.factor))

my_cforest_control = cforest_control(teststat = "quad", testtype = "Univ", mincriterion = 0, ntree = 500, mtry = 3, replace = FALSE)

m = cforest(as.factor(#{target_feature}) ~ ., data = data, control = my_cforest_control)
    

save(m, file = '#{self.file('model.Rdata')}')
data = NULL
    EOF

    confusion = data.R <<-EOF
load(file = '#{self.file('model.Rdata')}')
data = table(data$#{target_feature}, predict(m, OOB = TRUE))
    EOF

    Open.write(file('confusion.tsv'), confusion.to_s)
    confusion

#    importance = data.R <<-EOF
#rbbt.require('randomForest')
#load(file = '#{self.file('model.Rdata')}')
#data = m$importance
#    EOF
#
#    Open.write(file('importance.tsv'), importance.to_s)
#
#    err = data.R <<-EOF
#rbbt.require('randomForest')
#load(file = '#{self.file('model.Rdata')}')
#data = m$err.rate
#    EOF
#
#    log :err_rate, err[err.keys.last].first
#    Open.write(file('err_rate.tsv'), err.to_s)
#
#    importance
  end

  input :target_feature, :string, "name of feature to predict", 'anyCells'
  input :feature_files, :array, "Feature files to load", nil, :nofile => true
  task :svm_importance => :tsv do |target_feature,feature_files|
    data = Rbbt.data["Classification_influential_nodes.txt"].tsv(:fields => [target_feature], :type => :list)
    feature_files = Rbbt.data.glob("*.txt").collect{|f| f.find} if feature_files.nil?
    feature_files = feature_files - [Rbbt.data["Classification_influential_nodes.txt"].find]

    feature_files.each do |file|
      ftsv = TSV.open(file)
      ffields = ftsv.fields - ["HGNC.symbols"]
      ffields.reject!{|f| f != f.strip}
      data = data.attach ftsv, :fields => ffields
    end

    require 'rbbt/util/R'
    Open.mkdir files_dir
    data.R <<-EOF
rbbt.require('e1071')
names(data) = make.names(names(data))

data[is.na(data)] = -1
data = data[data$#{target_feature} != -1,]
m <- svm(as.factor(#{target_feature}) ~ ., data, na.action=na.omit)
save(m, file = '#{self.file('model.Rdata')}')
data = NULL
    EOF

    confusion = data.R <<-EOF
rbbt.require('randomForest')
load(file = '#{self.file('model.Rdata')}')
data = table(data$#{target_feature}, predict(m, OOB=T))
    EOF
    confusion
  end


  dep :importance, :compute => [:bootstrap, 5] do |jobname, options|
    jobs = []
    [true, false].each do |balance_data|
      %w(AGS COLO205 SW620 DU145 allCells anyCells).each do |target_feature|
        feature_files = Rbbt.data.glob("*.txt").collect{|f| f.find} if feature_files.nil?
        feature_files = feature_files - [Rbbt.data["Classification_influential_nodes.txt"].find]
        feature_files.each do |feature_file|
          jobname = [target_feature, (balance_data ? "balance_data" : "non_balanced_data"), File.basename(feature_file)] * "-"
          jobs << {:jobname => jobname, :inputs => options.merge(:balance_data => balance_data, :target_feature => target_feature, :feature_files => [feature_file])}
        end
        jobname = [target_feature, (balance_data ? "balance_data" : "non_balanced_data"), "all files"] * "-"
        jobs << {:jobname => jobname, :inputs => options.merge(:balance_data => balance_data, :target_feature => target_feature, :feature_files => nil)}
      end
    end
    jobs
  end
  task :batch => :tsv do
    importance = TSV.setup({}, "ID~#:type=:list")
    dependencies.each do |dep| 
      itsv = dep.file('importance.tsv').tsv
      itsv.fields = [dep.clean_name]
      importance = importance.attach itsv, :complete => true
    end

    Open.write(file('importance.tsv'), importance.to_s)

    tsv = TSV.setup({}, "Run~Target,Features,Balanced,Error#:type=:list")
    dependencies.each do |dep| 
      err_rate = dep.info[:err_rate]
      target_feature = dep.inputs[:target_feature]
      balance_data = dep.inputs[:balance_data]
      feature_files = (dep.inputs[:feature_files] || []).collect{|f| File.basename(f)}
      run = dep.clean_name
      tsv[run] = [target_feature, feature_files*",", balance_data, err_rate]
    end

    tsv
  end

end

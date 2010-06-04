#TODO this doesn't work yet - consider doing this later, use the install 
# -f|--force option for now

# distutils doesn't handle uninstalling things, this class deletes all the files
# this package installs if it has appropriate permissions to do it, otherwise
# print out the files that must be deleted to uninstall
class uninstall(build_py) :
  def run(self) :


    # delete modules
    print self.distribution.py_modules

    # delete extensions
    print self.distribution.ext_modules

    # delete packages
    print self.distribution.packages

    # delete package data
    print self.distribution.package_data

    # delete scripts
    print self.distribution.scripts

    print self.distribution.get_command_obj('install').get_outputs()

  def remove_path(self,path) :
    '''Attempt to remove the specified path, returning non-zero status code on error'''



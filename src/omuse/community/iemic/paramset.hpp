#ifndef PARAMSET_HPP
#define PARAMSET_HPP

#include <cstdint>
#include <typeindex>
#include <Teuchos_ParameterList.hpp>

class ParamSet
{
  public:
    template<typename T>
    struct tag {};

    template<class T>
    ParamSet(std::string name, tag<T>)
    : parameters(name)
    , defaultInitParams(T::getDefaultInitParameters())
    { parameters = defaultInitParams; }

    Teuchos::ParameterList& get();
    void reset();

    int get_num_params(const std::string&);
    std::string get_param_name(const std::string& param_name, int i);
    std::string get_param_type(const std::string& param_name);

    template<typename T>
    void
    get_param_value(const std::string& param_name, T& result)
    {
        const auto& param = get_param_entry(param_name);

        if (!param.isType<T>()) throw_param_mismatch(param_name);

        result = param.getValue(&result);
    }

    template<typename T>
    void
    set_param_value(const std::string& param_name, T result)
    {
        auto& param = get_param_entry(param_name);

        if (!param.isType<T>()) throw_param_mismatch(param_name);

        param.setValue(result);
    }

    template<typename T>
    void
    get_default_param_value(const std::string& param_name, T& result)
    {
        const auto& param = get_param_entry(param_name, defaultInitParams);

        if (!param.isType<T>()) throw_param_mismatch(param_name);

        result = param.getValue(& result);
    }

  private:
    Teuchos::ParameterEntry&
    get_param_entry(const std::string& param_name);

    static Teuchos::ParameterEntry&
    get_param_entry
    (const std::string& param_name, Teuchos::ParameterList& list);

    static const Teuchos::ParameterEntry&
    get_param_entry
    (const std::string& param_name, const Teuchos::ParameterList& list);

    static void throw_param_mismatch(const std::string& param_name);

    Teuchos::ParameterList parameters;
    const Teuchos::ParameterList defaultInitParams;

    static const std::map<std::type_index, std::string> types;
};
#endif

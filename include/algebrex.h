#ifndef ALGEBREX_H
#define ALGEBREX_H

#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <utility>
#include <vector>
#include <stack>

struct fraction_t {
	int up{};
	int down{1};

	fraction_t() = default;
	fraction_t(const int u, const int d) : up{u}, down{d} {
		for (int k = 2; k <= std::min(std::abs(up), std::abs(down)); ++k){
			while ((std::abs(up) % k == 0) && (std::abs(down) % k == 0)) {
				up /= k;
				down /= k;
			}
		}
	}
};

class value_t {
public:
	value_t() : fraction_value_(fraction_t{}) {}

	value_t(const fraction_t fract) {
		if (fract.down == 0) throw std::invalid_argument{"arithmatic error"};
		is_decimal_ = false;
		fraction_value_ = fract;
		decimal_value_ = static_cast<double>(fract.up) / fract.down;
		calculatable_ = true;
	}

	value_t(const double dv) {
		const auto new_value = value_t(std::to_string(dv));
		*this = new_value;
	}

	value_t(const std::string& str) {
		value_t newValue;
		if (const auto found = str.find('.'); found != std::string::npos) {
			const std::string left = str.substr(0, found);
			std::string right = str.substr(found + 1);

			if (left.length() >= 10)
				throw std::invalid_argument{"arithmatic error"};

			while (!right.empty() && right[right.length() -1] == '0')
				right = right.substr(0, right.length() -1);

			if (right.empty()) {
				newValue = value_t(fraction_t(stoi(left), 1));
				*this = newValue;
				return;
			}
			if (right.length() <= 5){
				if (right.find('.') != std::string::npos) throw std::invalid_argument{"arithmatic error"};
				const int leftNumber = stoi(left);
				const int rightNumber = stoi(right);
				const int multiplier = std::pow(10.0,(std::floor)(std::log10(rightNumber)) + 1.0);
				newValue = value_t(fraction_t(leftNumber*multiplier+rightNumber,multiplier));
				*this = newValue;
				return;
			}

			is_decimal_ = true;
			fraction_value_ = fraction_t();
			decimal_value_ = stod(str);
			calculatable_ = true;
			return;
		}
		newValue = value_t(fraction_t(stoi(str),1));
		*this = newValue;
	}

	[[nodiscard]] bool get_decimal() const { return is_decimal_; }
	[[nodiscard]] fraction_t getFracValue() const { return fraction_value_; }
	[[nodiscard]] double get_decimal_value() const { return decimal_value_; }
	[[nodiscard]] bool get_calculatable() const { return calculatable_; }

	// Print the value of the object
	[[nodiscard]]
	std::string print_value() const {
		if (!calculatable_) return {};
		if (is_decimal_) return std::to_string(decimal_value_);
		if (fraction_value_.down == 1 || fraction_value_.up == 0) return std::to_string(fraction_value_.up);

		return std::to_string(fraction_value_.up) + "/" + std::to_string(fraction_value_.down);
	}

	// Operator overload
	value_t& operator+=(const value_t value) {
		if (!calculatable_ || !value.get_calculatable()) { *this = value_t{}; return *this; }
		if (is_decimal_ || value.get_decimal()) *this = value_t{decimal_value_ + value.get_decimal_value()};
		else *this = value_t{fraction_t{
			fraction_value_.up*value.getFracValue().down + fraction_value_.down*value.getFracValue().up,
			value.getFracValue().down*fraction_value_.down}};
		return *this;
	}

	value_t& operator-=(const value_t value) {
		if (!calculatable_ || !value.get_calculatable()) { *this = value_t{}; return *this; }
		if (is_decimal_ || value.get_decimal()) *this = value_t(decimal_value_-value.get_decimal_value());
		else *this = value_t{fraction_t{
			fraction_value_.up*value.getFracValue().down - fraction_value_.down*value.getFracValue().up,
			value.getFracValue().down*fraction_value_.down}};
		return *this;
	}

	value_t& operator*=(const value_t z) {
		if (!calculatable_ || !z.get_calculatable()) { *this = value_t{}; return *this; }
		if (is_decimal_ || z.get_decimal()) *this = value_t{decimal_value_*z.get_decimal_value()};
		else *this = value_t{fraction_t{
			z.getFracValue().up * fraction_value_.up,
			z.getFracValue().down * fraction_value_.down}};
		return *this;
	}

	value_t& operator/=(const value_t z) {
		if (!calculatable_ || !z.get_calculatable()) { *this = value_t{}; return *this; }
		if (is_decimal_ || z.get_decimal()) *this = value_t(decimal_value_/z.get_decimal_value());
		else *this = value_t{fraction_t{
			z.getFracValue().down * fraction_value_.up,
			z.getFracValue().up*fraction_value_.down}};
		return *this;
	}

	value_t& powv(const value_t value) {
		if (!calculatable_ || !value.get_calculatable()) { *this = value_t(); return *this; }

		if ((decimal_value_ < 0) && (value.get_decimal_value() != static_cast<double>(static_cast<int>(value.get_decimal_value()))))
			throw std::invalid_argument{"arithmatic error"};

		if (!is_decimal_ && value.get_decimal_value() == static_cast<double>(static_cast<int>(value.get_decimal_value()))) {
			fraction_value_.up = std::pow(fraction_value_.up, fabs(value.get_decimal_value()));
			fraction_value_.down = std::pow(fraction_value_.down, fabs(value.get_decimal_value()));

			if (value.get_decimal_value() < 0) {
				std::swap(fraction_value_.up,fraction_value_.down);
				if (fraction_value_.down < 0) {
					fraction_value_.up *= -1;
					fraction_value_.down *= -1;
				}
			}
		} else { *this = value_t(pow(decimal_value_,value.get_decimal_value())); }
		return *this;
	}

	friend value_t operator+(value_t a, const value_t b) { return a += b; }
	friend value_t operator-(value_t a, const value_t b) { return a -= b; }
	friend value_t operator*(value_t a, const value_t b) { return a *= b; }
	friend value_t operator/(value_t a, const value_t b) { return a /= b; }
	friend value_t operator-(value_t a) { return a *= -1; }
	friend value_t powv(value_t a, value_t b) { return a.powv(b); }
	friend std::ostream &operator<<(std::ostream &out, const value_t &m) {
		return out << m.print_value();
	}
private:
	bool is_decimal_{};
	fraction_t fraction_value_;
	double decimal_value_{};
	bool calculatable_{};
};

struct variable_t {
	std::string name;
	value_t value;

	variable_t(std::string nm, const value_t &val) : name{std::move(nm)}, value{val} {}
};

struct function_t {
	std::string name;
	double (*func)(double);

	function_t(std::string nm, double (*f)(double)) : name{std::move(nm)}, func{f} {}
};

enum class block_type { number, symbol, func, constant, var, opening_bracket, closing_bracket, null };

struct block_t {
	int start{};
	int end{};
	int level{};
	block_type type{block_type::null};

	block_t() = default;
	block_t(const int s, const int e, const int lvl, const block_type ty)
		: start{s}, end{e}, level{lvl}, type{ty} {}
};

class expr_solver {
public:
	expr_solver() { add_predefined(); }

	std::string solve(std::string expr) {
		expr = discard_spaces(expr);
		bool is_declaration{};
		std::string new_var_name;

		if (expr.empty()) {
			throw std::invalid_argument{"invalid expression"};
		}
		const bool declaration_valid = check_declaration(expr, new_var_name, is_declaration);
		if (!declaration_valid) {
			blocks_.clear();
			throw std::invalid_argument{"declaration invalid"};
		}
		deal_with_negative_sign(expr);

		if (const bool group_succeed = group_expr(expr); !group_succeed) {
			blocks_.clear();
			throw std::invalid_argument{"group not succeed"};
		}
		const value_t result = calculate_expr(expr, 0, blocks_.size());
		std::string output_string;

		if (result.get_calculatable()) {
			if (is_declaration) {
				for (auto& constant : constants_) {
					if (new_var_name == constant.name) {
						blocks_.clear();
						throw std::logic_error{"constant " + new_var_name + "cannot be declared"};
					}
				}
				bool variable_declared_before{};

				for (auto& variable : variables_) {
					if (new_var_name == variable.name) {
						variable_declared_before = true;
						variable = variable_t(new_var_name, result);
					}
				}
				if (!variable_declared_before) {
					variables_.emplace_back(new_var_name, result);
				}
				std::cout << new_var_name << " = " << result.print_value();
			} else {
				constants_[constants_.size()-1] = variable_t("ans", result);
				output_string = result.print_value();
			}
		} else {
			blocks_.clear();
			throw std::logic_error{"error"};
		}
		blocks_.clear();
		return output_string;
	}

private:
	void add_predefined() {
		constants_.emplace_back("pi", value_t(M_PI));
		constants_.emplace_back("ans", value_t());
		functions_.emplace_back("exp", exp);
		functions_.emplace_back("sqrt",sqrt);
		functions_.emplace_back("floor", floor);
		functions_.emplace_back("ln", log);
		functions_.emplace_back("log", log10);
	}

	static std::string discard_spaces(const std::string &str) {
		std::string new_str;
		for (std::size_t i{}; i < str.length(); ++i)
			if (!isspace(str[i])) new_str += str[i];
		return new_str;
	}

	static bool check_declaration(std::string &expr, std::string &new_var_name, bool &is_dec) {
		if (const auto found = expr.find('='); found != std::string::npos) {
			new_var_name = expr.substr(0,found);
			is_dec = true;
			expr = expr.substr(found + 1);

			if (expr.find('=') != std::string::npos) return false;
			if (new_var_name.empty() || !isalpha(new_var_name[0])) return false;

			for (int i = 1; i < new_var_name.length(); i++)
				if (!(isalnum(new_var_name[i]) || new_var_name[i] == '_')) return false;
		}
		return true;
	}

	bool group_expr(const std::string &expr) {
		auto start{0};
		auto level{0};
		auto current_type = block_type::null;

		for (int i = 0; i <= expr.length(); ++i) {
			const block_type this_type = char_type(expr[i]);
			bool need_new_block{};
			need_new_block |= (current_type == block_type::opening_bracket);
			need_new_block |= (current_type == block_type::closing_bracket);
			need_new_block |= (current_type == block_type::symbol);
			need_new_block |= (this_type != current_type);
			need_new_block |= (i == expr.length());
			need_new_block &= !(this_type == block_type::number && current_type == block_type::func);

			if (need_new_block) {
				if (current_type == block_type::func) {
					current_type = analyze_str_type(expr.substr(start, i-start));
					if (current_type == block_type::null) {
						return false;
					}
				}
				if (i != 0) {
					auto new_block = block_t{start, i, level, current_type};
					blocks_.push_back(new_block);
				}
				if (current_type == block_type::closing_bracket) {
					--level;
				}
				current_type = this_type;
				start = i;
			}
			if (this_type == block_type::opening_bracket) {
				++level;
			}
		}
		if (level != 0) {
			return false;
		}
		return true;
	}

	[[nodiscard]]
	block_type analyze_str_type(const std::string &str) const {
		for (const auto & function : functions_)
			if (str == function.name) return block_type::func;

		for (const auto & i : constants_)
			if (str == i.name) return block_type::constant;

		for (const auto & variable : variables_)
			if (str == variable.name) return block_type::var;

		return block_type::null;
	}

	static block_type char_type(const char c) {
		if (c == '_' || isalpha(c)) return block_type::func;
		if (c == '.' || isdigit(c)) return block_type::number;
		if (c == '(') return block_type::opening_bracket;
		if (c == ')') return block_type::closing_bracket;
		if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') return block_type::symbol;
		return block_type::null;
	}

	static void deal_with_negative_sign(std::string &exp) {
		if (!exp.empty() && exp[0] == '-') exp = '0' + exp;

		for (int i = 1; i < exp.length(); i++)
			if (exp[i] == '-' && exp[i - 1] == '(') exp = exp.substr(0, i) + '0' + exp.substr(i);
	}

	value_t calculate_expr(const std::string& exp, int start_block, int end_block) {
		std::stack<value_t> values;
		std::stack<char> ops;

		for (int i = end_block-1; i >= start_block; i--) {
			std::string block_str = exp.substr(blocks_[i].start, blocks_[i].end-blocks_[i].start);
			int i_increment{};

			if (blocks_[i].type == block_type::number) {
				values.emplace(block_str);
			} else if (blocks_[i].type == block_type::func) {
				return value_t{};
			} else if (blocks_[i].type == block_type::constant) {
				value_t const_value;
				for (auto& constant : constants_) {
					if (block_str == constant.name) {
						const_value = constant.value;
						if (!const_value.get_calculatable()) return {};
						break;
					}
				}
				values.push(const_value);
			} else if (blocks_[i].type == block_type::var) {
				value_t var_value;
				for (auto & variable : variables_) {
					if (block_str == variable.name) {
						var_value = variable.value;
						break;
					}
				}
				values.push(var_value);
			} else if (blocks_[i].type == block_type::closing_bracket) {
				if (auto cor_block = find_index_of_bracket_ending(i);
					(cor_block != 0) && (blocks_[cor_block - 1].type == block_type::func))
				{
					double (*func_to_use)(double);
					std::string func_name = exp.substr(
						blocks_[cor_block-1].start,
						blocks_[cor_block-1].end - blocks_[cor_block - 1].start);

					for (auto & function : functions_)
						if (func_name == function.name) { func_to_use = function.func; break; }

					value_t value_in_func = calculate_expr(exp, cor_block+1, i);

					if (!value_in_func.get_calculatable()) return {};
					if (value_in_func.get_decimal_value() < 0 && func_name == "sqrt") throw std::invalid_argument{"arithmatic error"};

					double func_result = (*func_to_use)(value_in_func.get_decimal_value());
					auto new_value = value_t(func_result);
					values.push(new_value);

					i_increment -= i - cor_block + 1;
				} else {
					value_t value_to_push = calculate_expr(exp, cor_block + 1, i);
					if (!value_to_push.get_calculatable()) return {};

					values.push(value_to_push);
					i_increment -= i - cor_block;
				}
			} else if (blocks_[i].type == block_type::symbol) {
				if (block_str[0] == '*' || block_str[0] == '/') {
					if (values.size() != ops.size() + 1) return {};

					while(!ops.empty()) {
						if (char last_op = ops.top(); last_op == '^') {
							ops.pop();
							value_t op1 = values.top(); values.pop();
							value_t op2 = values.top(); values.pop();
							values.push(powv(op1,op2));
						} else break;
					}
				} else if ((block_str[0] == '+') || (block_str[0] == '-')) {
					if (values.size() != ops.size() + 1) return {};

					while(!ops.empty()) {
						if (char last_op = ops.top(); !(last_op == '+' || last_op == '-')) {
							ops.pop();
							value_t op1 = values.top(); values.pop();
							value_t op2 = values.top(); values.pop();

							if (last_op == '*') {
								values.push(op1*op2);
							} else if(last_op == '/') {
								values.push(op1/op2);
							} else if(last_op == '^') {
								values.push(powv(op1,op2));
							}
						} else break;
					}
				}
				ops.push(block_str[0]);
			} else return {};
			i += i_increment;
		}
		if (values.size() != (ops.size() + 1)) throw std::invalid_argument{"invalid expression"};
		while(!ops.empty()) {
			char last_op = ops.top();
			ops.pop();
			value_t op1 = values.top(); values.pop();
			value_t op2 = values.top(); values.pop();

			if (last_op == '+') values.push(op1 + op2);
			else if (last_op == '-') values.push(op1-op2);
			else if (last_op == '*') values.push(op1*op2);
			else if (last_op == '/') values.push(op1/op2);
			else if (last_op == '^') values.push(powv(op1, op2));
		}
		if (values.size() > 1) throw std::invalid_argument{"invalid expression"};
		value_t return_value = values.top();

		while(!ops.empty()) ops.pop();
		while(!values.empty()) values.pop();

		return return_value;
	}

	[[nodiscard]]
	int find_index_of_bracket_ending(const int block_id) const {
		const auto level_to_find = blocks_[block_id].level - 1;
		auto current_block_id = block_id;

		while (blocks_[current_block_id].level != level_to_find) {
			if (current_block_id == 0) return current_block_id;
			current_block_id--;
		}
		return current_block_id + 1;
	}

private:
	std::vector<block_t> blocks_;
	std::vector<variable_t> variables_;
	std::vector<variable_t> constants_;
	std::vector<function_t> functions_;
};

#endif //ALGEBREX_H

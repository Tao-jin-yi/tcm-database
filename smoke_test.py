"""新网站冒烟测试：AppTest 真机运行 app.py，检查整页异常。"""
import sys

from streamlit.testing.v1 import AppTest


def main() -> int:
    at = AppTest.from_file("app.py", default_timeout=300)
    at.run()
    if at.exception:
        print("异常:", at.exception[0].value)
        return 1
    print("页面加载 OK")

    # 交互：中药搜索
    at.text_input[0].set_value("甘草").run()
    if at.exception:
        print("中药搜索异常:", at.exception[0].value)
        return 1
    print("中药搜索 OK")

    # 交互：成分搜索
    at.text_input[1].set_value("quercetin").run()
    if at.exception:
        print("成分搜索异常:", at.exception[0].value)
        return 1
    print("成分搜索 OK")

    # 交互：疾病搜索
    at.text_input[3].set_value("diabetes").run()
    if at.exception:
        print("疾病搜索异常:", at.exception[0].value)
        return 1
    print("疾病搜索 OK")
    print("全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())

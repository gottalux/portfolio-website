<template>
  <v-app>
    <v-navigation-drawer app right permanent color="secondary" dark>
      <template v-slot:prepend>
        <v-list-item two-line>
          <v-list-item-avatar size="224">
            <img src="./assets/avatar.jpg" />
          </v-list-item-avatar>
        </v-list-item>
        <v-list-item three-line color="primary" dark>
          <v-list-item-content>
            <v-list-item-title>Adahlia S. Lux</v-list-item-title>
            <v-list-item-subtitle class="primary--text"
              >Final year PhD student, supervised by
              <a href="http://www.cs.bham.ac.uk/~axk/" target="_blank"
                >Prof. Ata Kab&aacute;n
              </a></v-list-item-subtitle
            >
          </v-list-item-content>
        </v-list-item>
      </template>
      <v-divider color="accent"></v-divider>

      <v-list dense>
        <template v-for="(item, index) in navigationItems">
          <v-divider
            v-if="item.divider"
            :key="index"
            color="accent"
            inset
          ></v-divider>
          <v-list-item v-else :key="index" color="primary" dark>
            <v-list-item-icon>
              <v-icon>{{ item.icon }}</v-icon>
            </v-list-item-icon>

            <v-list-item-content>
              <v-list-item-title>{{ item.title }}</v-list-item-title>
            </v-list-item-content>
          </v-list-item>
        </template>
      </v-list>

      <template v-slot:append>
        <v-divider color="accent"></v-divider>

        <v-list dense>
          <template v-for="(item, index) in externalLinkItems">
            <v-divider
              v-if="item.divider"
              :key="index"
              color="accent"
              inset
            ></v-divider>
            <v-list-item
              v-else
              :key="index"
              color="primary"
              dark
              @click="goToUrl(item.url)"
            >
              <v-list-item-icon>
                <img
                  v-if="item.icon === 'linkedin'"
                  src="./assets/linkedin.png"
                />
                <v-icon v-else>{{ item.icon }}</v-icon>
              </v-list-item-icon>

              <v-list-item-content>
                <v-list-item-title>{{ item.title }}</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
          </template>
        </v-list>
      </template>
    </v-navigation-drawer>

    <v-main>
      <HelloWorld />
    </v-main>
  </v-app>
</template>

<script>
import HelloWorld from "./components/HelloWorld";

export default {
  name: "App",

  components: {
    HelloWorld,
  },

  metaInfo: {
    title: "Home Page",
    titleTemplate: "%s - Adahlia Lux",
  },
  data: () => ({
    navigationItems: [
      {
        icon: "engineering",
        title: "Research interests",
      },
      { divider: true },
      {
        icon: "newspaper",
        title: "Publications",
      },
      { divider: true },
      {
        icon: "groups",
        title: "Demonstrating",
      },
      { divider: true },
      {
        icon: "co_present",
        title: "Slides",
      },
      { divider: true },
      {
        icon: "portrait",
        title: "Background",
      },
      { divider: true },
      {
        icon: "outgoing_mail",
        title: "Contact me",
      },
    ],
    externalLinkItems: [
      {
        icon: "palette",
        title: "Personal and Art",
        url: "https://allmylinks.com/gottalux",
      },
      { divider: true },
      {
        icon: "linkedin",
        title: "LinkedIn",
        url: "https://www.linkedin.com/in/adahlia-lux-27121075/",
      },
    ],
  }),
  methods: {
    goToUrl(url) {
      window.open(url, "_blank");
    },
  },
};
</script>
